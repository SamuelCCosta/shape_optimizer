from abc import ABC, abstractmethod
import numpy as np
import random
import multiprocessing as mp
from pebble import concurrent, ProcessExpired

spawn_ctx = mp.get_context('spawn') #fork context (default on linux) can lead to some issues

class BaseSimmulatedAnnealing(ABC):
    '''Regular Simmulated Annealing algorithm, with geometric temperature reduction'''

    def __init__(self, initial_temp = 1000.0, min_temp = 0.001, cooling_rate = 0.95):
        self.initial_temp = initial_temp
        self.temp = initial_temp
        self.min_temp = min_temp
        self.cooling_rate = cooling_rate
        #track best state and energy
        self.best_state = None
        self.best_cost = float('inf')

    @abstractmethod
    def get_neighbour(self, state):
        '''Given a state, get a neighbouring state.'''
        pass

    @abstractmethod
    def raw_cost_function(self,state) -> float:
        '''Get cost, might be infinite, error out or timeout'''
        pass

    @concurrent.process(timeout=2.0, mp_context=spawn_ctx)
    def _async_cost(self, state):
        '''async wrapper for the unstable raw_cost_function'''
        return self.raw_cost_function(state)

    def get_cost(self, state) -> float:
        '''Safely return the cost functional, or infinite if it did not compute'''
        try:
            future = self._async_cost(state)
            return future.result() #type: ignore #Pylance is not smart enough to know it returns a Future object
        except:
            return float('inf')
    
    def acceptance_probability(self, old_cost, new_cost):
        '''Calculate the probability of accepting the solution.'''
        if new_cost < old_cost:
            return 1.0
        return np.exp((old_cost - new_cost) / self.temp)
    
    def update_temperature(self):
        '''Updates temperature. Defaults to fixed geometric cooling.'''
        self.temp *= self.cooling_rate
        return None
    
    def run(self, initial_state):
        '''The main optimization loop, using the rules defined in the Class'''
        current_state = initial_state
        current_cost = self.get_cost(initial_state)
        assert(current_cost != float('inf'))

        self.best_state = current_state
        self.best_cost = current_cost

        while self.temp > self.min_temp:
            neighbour_cost = float('inf')
            neighbour_state = None
            while neighbour_cost == float('inf'):
                neighbour_state = self.get_neighbour(current_state)
                neighbour_cost = self.get_cost(neighbour_state)
            #here we have suitable parameters and cost

            if random.random() < self.acceptance_probability(current_cost, neighbour_cost):
                current_state = neighbour_state
                current_cost = neighbour_cost

                if current_cost < self.best_cost:
                    self.best_state = current_state
                    self.best_cost = current_cost
            
            self.update_temperature()
        
        return self.best_state, self.best_cost


if __name__ == '__main__':
    #example
    class ParabolaMinimizer(BaseSimmulatedAnnealing):
        def __init__(self, step_size = 1.0, **kwargs):
            super().__init__(**kwargs)
            self.step_size = step_size

        def raw_cost_function(self, state) -> float:
            return state ** 2
        
        def get_neighbour(self, state):
            size = self.step_size * self.temp / self.initial_temp
            noise = random.normalvariate(0, size)
            return state + noise
        
    solver = ParabolaMinimizer(step_size = 50.0,
                               initial_temp = 500.0,
                               min_temp = 1e-04,
                               cooling_rate = 0.99)
    
    initial_guess = 100.0
    best_param, lowest_cost = solver.run(initial_guess)

    print(f'{best_param= :.4f},{lowest_cost= :.2e}')
