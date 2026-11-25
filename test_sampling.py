import copy
import json
import random
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from db_utils import get_column_names_dict
from orchestrator import (
    EllipseSA, EllipseDSA, build_grid_combinations, build_zip_combinations,
    build_lhs_combinations, run_experiments,
)
from shape_optimizer import shape_optimizer


GEOMETRY = {
    'geometric_config': {'x_max': 1.0, 'y_max': 1.0},
    'h': 0.02, 'num_ellipses': 2,
}
PENALIZATIONS = {'linear': 1.0}
SCALES = (0.1, 0.2, 1.0, 2.0, 0.0)


def sampling_worker(queue, run_id, geometry, penalizations, kwargs):
    solver = EllipseSA(geometry, penalizations, **kwargs)
    queue.put({'run_id': run_id, 'best_cost': 0.0,
               'best_param': json.dumps(solver.get_neighbour([0.0] * 10))})


class SamplingTests(unittest.TestCase):
    def make_solver(self, solver_type, **kwargs):
        key = 'perturbation' if solver_type is EllipseSA else 'initial_perturbation'
        return solver_type(GEOMETRY, PENALIZATIONS, **{key: SCALES}, **kwargs)

    def test_default_and_invalid_sampling(self):
        for solver_type in (EllipseSA, EllipseDSA):
            with self.subTest(solver=solver_type.__name__):
                self.assertEqual(self.make_solver(solver_type).sampling, 'normal')
                for invalid in ('gaussian', '', None, ['normal']):
                    with self.assertRaisesRegex(ValueError, 'sampling'):
                        self.make_solver(solver_type, sampling=invalid)

    def test_step_distributions_and_scales(self):
        original_rng = random.getstate()
        self.addCleanup(random.setstate, original_rng)
        for solver_type in (EllipseSA, EllipseDSA):
            for sampling, expected_fraction in (('normal', 0.6827), ('uniform', 1.0), ('cauchy', 0.5)):
                with self.subTest(solver=solver_type.__name__, sampling=sampling):
                    random.seed(127)
                    solver = self.make_solver(solver_type, sampling=sampling)
                    state = [10.0] * 10
                    steps = [solver.get_neighbour(state) for _ in range(10000)]
                    for index, scale in enumerate(SCALES * 2):
                        if scale == 0.0:
                            self.assertTrue(all(row[index] == 10.0 for row in steps))
                            continue
                        normalized = [(row[index] - 10.0) / scale for row in steps]
                        fraction = sum(abs(step) <= 1.0 + 1e-12 for step in normalized) / len(steps)
                        self.assertAlmostEqual(fraction, expected_fraction, delta=0.025)
                        positive = sum(step > 0 for step in normalized) / len(steps)
                        self.assertAlmostEqual(positive, 0.5, delta=0.025)
                        if sampling == 'cauchy':
                            tail = sum(abs(step) > 4 for step in normalized) / len(steps)
                            self.assertAlmostEqual(tail, 0.156, delta=0.025)
                    self.assertEqual(state, [10.0] * 10)

    def test_dsa_uses_updated_perturbation(self):
        solver = self.make_solver(EllipseDSA, sampling='uniform', perturbation_update=0.5)
        solver.update_perturbation()
        for _ in range(100):
            neighbour = solver.get_neighbour([0.0] * 10)
            self.assertTrue(all(abs(step) <= scale / 2 for step, scale in zip(neighbour, SCALES * 2)))

    def test_sampling_sweeps(self):
        kwargs = {'perturbation': SCALES, 'sampling': ['normal', 'uniform', 'cauchy']}
        for builder in (build_grid_combinations, build_zip_combinations, build_lhs_combinations):
            with self.subTest(builder=builder.__name__):
                extra = {'num_samples': 3} if builder is build_lhs_combinations else {}
                bundles = builder(GEOMETRY, PENALIZATIONS, kwargs, **extra)
                self.assertCountEqual([b['kwargs_optimization']['sampling'] for b in bundles], kwargs['sampling'])

    def test_sql_column_order_independent_of_config_order(self):
        for perturbation in ('perturbation', 'initial_perturbation'):
            columns = get_column_names_dict(GEOMETRY, PENALIZATIONS, {
                'sampling': 'uniform', perturbation: SCALES, 'gen_seed': 12,
            })
            names = list(columns)
            index = names.index('gen_seed')
            self.assertEqual(names[index:index + 3], ['gen_seed', 'sampling', perturbation])
            self.assertEqual(columns['sampling'], 'TEXT')

    def test_json_config_sampling_reaches_experiments(self):
        config = {
            'optimization_type': 'SA', 'combination_type': 'grid', 'n_runs': 1,
            'geometric_params': GEOMETRY, 'penalizations': PENALIZATIONS,
            'kwargs_SA': {'perturbation': SCALES, 'sampling': ['normal', 'uniform', 'cauchy']},
            'database_name': 'unused.db', 'table_name': 'results', 'max_processes': 1,
        }
        for optimization in ('SA', 'DSA'):
            with self.subTest(optimization=optimization), tempfile.TemporaryDirectory() as directory:
                variant = copy.deepcopy(config)
                variant['optimization_type'] = optimization
                if optimization == 'DSA':
                    variant['kwargs_DSA'] = variant.pop('kwargs_SA')
                    variant['kwargs_DSA']['initial_perturbation'] = variant['kwargs_DSA'].pop('perturbation')
                path = Path(directory) / 'config.json'
                path.write_text(json.dumps(variant))
                with patch('shape_optimizer.run_experiments') as run:
                    shape_optimizer(str(path))
                bundles = run.call_args.args[1]
                self.assertEqual([b['kwargs_optimization']['sampling'] for b in bundles], ['normal', 'uniform', 'cauchy'])

    def test_experiments_persist_selected_and_default_sampling(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / 'runs.db')
            kwargs = {'perturbation': SCALES, 'sampling': ['normal', 'uniform', 'cauchy']}
            bundles = build_grid_combinations(GEOMETRY, PENALIZATIONS, kwargs)
            run_experiments(sampling_worker, bundles, db_path=path, max_processes=3)
            default_bundles = build_grid_combinations(GEOMETRY, PENALIZATIONS, {'perturbation': SCALES})
            run_experiments(sampling_worker, default_bundles, db_path=path)
            with sqlite3.connect(path) as conn:
                rows = conn.execute('SELECT sampling, best_cost, best_param FROM results ORDER BY run_id').fetchall()
                self.assertEqual([r[0] for r in rows], ['normal', 'uniform', 'cauchy', 'normal'])
                self.assertTrue(all(r[1] == 0.0 and len(json.loads(r[2])) == 10 for r in rows))
                names = [r[1] for r in conn.execute('PRAGMA table_info(results)')]
                index = names.index('gen_seed')
                self.assertEqual(names[index:index + 3], ['gen_seed', 'sampling', 'perturbation'])

    def test_invalid_sampling_fails_before_creating_database(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'invalid.db'
            bundles = build_grid_combinations(GEOMETRY, PENALIZATIONS, {
                'perturbation': SCALES, 'sampling': ['normal', 'invalid'],
            })
            with self.assertRaisesRegex(ValueError, 'sampling'):
                run_experiments(sampling_worker, bundles, db_path=str(path))
            self.assertFalse(path.exists())

    def test_legacy_table_is_unchanged_when_using_a_new_table(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / 'runs.db')
            with sqlite3.connect(path) as conn:
                conn.execute('CREATE TABLE results (run_id INTEGER PRIMARY KEY AUTOINCREMENT, perturbation TEXT, gen_seed INTEGER)')
                conn.execute("INSERT INTO results (perturbation, gen_seed) VALUES ('[0.1]', 7)")
            bundles = build_grid_combinations(GEOMETRY, PENALIZATIONS, {
                'perturbation': SCALES, 'sampling': 'uniform',
            })
            run_experiments(sampling_worker, bundles, db_path=path, table_name='new_results')
            with sqlite3.connect(path) as conn:
                names = [r[1] for r in conn.execute('PRAGMA table_info(results)')]
                self.assertEqual(names, ['run_id', 'perturbation', 'gen_seed'])
                self.assertEqual(conn.execute('SELECT * FROM results').fetchall(), [(1, '[0.1]', 7)])
                self.assertEqual(conn.execute('SELECT sampling FROM new_results').fetchall(), [('uniform',)])
                self.assertEqual(conn.execute('PRAGMA integrity_check').fetchone(), ('ok',))

    def test_targeting_a_legacy_table_does_not_alter_it(self):
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / 'runs.db')
            with sqlite3.connect(path) as conn:
                conn.execute('CREATE TABLE results (run_id INTEGER PRIMARY KEY AUTOINCREMENT, perturbation TEXT, gen_seed INTEGER)')
                conn.execute("INSERT INTO results (perturbation, gen_seed) VALUES ('[0.1]', 7)")
            bundles = build_grid_combinations(GEOMETRY, PENALIZATIONS, {
                'perturbation': SCALES, 'sampling': 'uniform',
            })
            with self.assertRaisesRegex(sqlite3.OperationalError, 'sampling'):
                run_experiments(sampling_worker, bundles, db_path=path)
            with sqlite3.connect(path) as conn:
                names = [r[1] for r in conn.execute('PRAGMA table_info(results)')]
                self.assertEqual(names, ['run_id', 'perturbation', 'gen_seed'])
                self.assertEqual(conn.execute('SELECT * FROM results').fetchall(), [(1, '[0.1]', 7)])
                self.assertEqual(conn.execute('PRAGMA integrity_check').fetchone(), ('ok',))


if __name__ == '__main__':
    unittest.main()
