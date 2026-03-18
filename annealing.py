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

class DirectSimmulatedAnnealing(ABC):
    def __init__(self, initial_temp = 1000.0, min_temp = 0.0001, cooling_rate_max = 0.999,
                 cooling_rate_min = 0.99, base_markov_length = 50, num_configs = 100,
                 initial_perturbation = 0.05, min_perturbation = 0.01, min_cost_gap = 0.0001):
        self.initial_temp = initial_temp
        self.temp = initial_temp
        self.min_temp = min_temp
        self.min_cost_gap = min_cost_gap

        #cooling
        self.max_cooling_rate = cooling_rate_max
        self.min_cooling_rate = cooling_rate_min
        self.effective_cooling = 1.0

        #markov lengths
        self.base_markov_length = base_markov_length
        self.current_markov_length = 0
        self.effective_past_markov_length = 0

        self.num_configs = num_configs

        #perturbation values
        self.perturbation = initial_perturbation
        self.min_perturbation = min_perturbation
        
        #track best state and energy
        self.configurations = [] #sorted by cost

    @abstractmethod
    def get_neighbour(self, state):
        '''Given a state, get a neighbouring state.'''
        pass

    @abstractmethod
    def raw_cost_function(self,state) -> float:
        '''Get cost, might be infinite, error out or timeout'''
        pass

    @abstractmethod
    def get_starting_configs(self):
        '''Get starting configuration to start the optimization process'''
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
    
    def update_temperature(self, effective_markov_length):
        '''Updates temperature, which depends on the length of past and current Markov Chain, as well as past cooling parameter.'''
        cooling_rate = 1.0
        if self.effective_past_markov_length == 0:
            cooling_rate = (self.min_cooling_rate + self.max_cooling_rate) / 2.0
        elif self.current_markov_length > effective_markov_length:
            cooling_rate = self.max_cooling_rate
        elif self.current_markov_length > self.effective_past_markov_length:
            cooling_rate = self.effective_cooling - (self.effective_cooling - self.min_cooling_rate) * (1 - self.effective_past_markov_length/self.current_markov_length)
        else:
            cooling_rate = self.max_cooling_rate - (self.max_cooling_rate - self.effective_cooling) * (self.current_markov_length/self.effective_past_markov_length)
        #update temperature
        self.effective_cooling = cooling_rate
        self.temp *= self.effective_cooling
        return None
    
    def update_perturbation(self):
        '''Updates how much perturbation is allowed. Default: geometric reduction.'''
        self.perturbation *= 0.995
        return None

    def cost_gap(self):
        '''Returns the difference between highest and lowest cost'''
        return self.configurations[-1][0] - self.configurations[0][0]

    def markov_length(self):
        '''Gives the maximum Markov chain length dependant in the gap in configurations cost.'''
        return self.base_markov_length * (2 - np.exp(-self.cost_gap()))

    def best_cost(self):
        '''Returns lowest cost'''
        return self.configurations[0][0]
    
    def worst_cost(self):
        return self.configurations[-1][0]
    
    def second_worst_cost(self):
        '''AKA worse cost'''
        return self.configurations[-2][0]
    
    def update_configuration(self, new_cost, new_state):
        '''Removes the worst configuration and adds the new_state, and then orders the list again.'''
        self.configurations.pop() #remove worst
        self.configurations.append((new_cost, new_state))
        self.configurations.sort() #could just insert in the right place
        return None

    def run(self):
        '''The main optimization loop, returns a tuple with best cost and best configuration.'''
        while (self.perturbation > self.min_perturbation) or (self.temp > self.min_temp) or (self.cost_gap() >  self.min_cost_gap):
            self.current_markov_length = self.markov_length()
            previous_second_worst_cost = self.second_worst_cost()
            length = 0
            while length < self.current_markov_length: 
                length += 1
                neighbour_state = None
                neighbour_cost = float('inf')
                
                #get suitable configuration
                while neighbour_cost == float('inf'):
                    chosen_cost, chosen_state = random.choice(self.configurations)
                    neighbour_state = self.get_neighbour(chosen_state)
                    neighbour_cost = self.get_cost(neighbour_state)
                
                if random.random() < self.acceptance_probability(self.worst_cost(),neighbour_cost):
                    current_best_cost = self.best_cost()
                    self.update_configuration(neighbour_cost, neighbour_state)

                    if neighbour_cost < current_best_cost:
                        break
                
            self.update_temperature(length)
            if self.second_worst_cost() >= previous_second_worst_cost: #worse configuration not improved
                self.update_perturbation()
            self.effective_past_markov_length = length
        
        return self.configurations[0]



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
                               min_temp = 0.001,
                               cooling_rate = 0.99)
    
    initial_guess = 100.0
    best_param, lowest_cost = solver.run(initial_guess)

    print(f'{best_param= :.4f},{lowest_cost= :.2e}')

