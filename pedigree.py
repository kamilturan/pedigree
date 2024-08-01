import copy
import random
from enum import Enum
from typing import (Any, Optional, Self, List)
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

class Gender(Enum):
    Female = 0
    Male = 1
    
class Status(Enum):
    Healthy = 0
    Affected = 1

class Somatic:
    
    def __init__(self) -> None:
        self.__genotype = ""

    @property
    def allele(self) -> str:
        return self.__genotype

    @allele.setter
    def allele(self, value) -> None:
        self.__genotype = value
        
    def choice(self):
        return random.choice(self.__genotype)
    
    def __repr__(self) -> str:
        return self.allele
    
class Person:
    
    def __init__(self, id: int,
                 gender: Gender,
                 status: Status,
                 father: Optional[Self] = None,
                 mother: Optional[Self] =None) -> None:
        self.__id = id
        self.__gender = gender
        self.__status = status
        self.__father = father
        self.__mother = mother
        self.__soma = Somatic()
    
    @property
    def somatic(self) -> Somatic:
        return self.__soma
    
    @somatic.setter
    def somatic(self, value) -> None:
        self.__soma = value
        
    @property 
    def id(self) -> int:
        return self.__id
    
    @property
    def gender(self) -> Gender:
        return self.__gender
    
    @property
    def status(self) -> Status:
        return self.__status
    
    @property
    def father(self) -> Self:
        return self.__father
    
    @property
    def mother(self) -> Self:
        return self.__mother
    
    def inherit(self):
        if (self.mother is not None) and (self.father is not None):
            assert self.mother.somatic.allele != "", "Mother somatic alleles are not none."
            assert self.father.somatic.allele != "", "Father somatic alleles are not none."
            a1 = self.mother.somatic.choice()
            a2 = self.father.somatic.choice()
            self.somatic.allele = f"{a1}{a2}"
    
    def __eq__(self, other: Self) -> bool:
        return (self.id == other.id) and (self.gender == other.gender) and (self.status == other.status)
    
    def __repr__(self) -> str:
        return f"[{self.id}] {self.gender.name} {self.status.name} S:{self.__soma}"

class Models:
    
    AutosomalDominant = {Status.Affected:["AA","Aa","aA"],
                         Status.Healthy:["aa"],
                         "code":"AD",
                         "name":"Autosomal dominant inheritance",
                         "reverse":{"AA":Status.Affected,
                                    "Aa":Status.Affected,
                                    "aA":Status.Affected,
                                    "aa":Status.Healthy} }
    AutosomalRecessive = {Status.Affected:["aa"],
                          Status.Healthy:["AA","Aa","aA"],
                          "code":"AR",
                          "name":"Autosomal dominant inheritance",
                          "reverse":{"AA":Status.Healthy,
                                    "Aa":Status.Healthy,
                                    "aA":Status.Healthy,
                                    "aa":Status.Affected} }
    
    YLinkedInheritance = {Status.Affected:["A","a"],
                          Status.Healthy:["-"],
                          "code":"YL",
                          "name":"Y linked inheritance",
                          "reverse":{"a":Status.Affected,
                                     "A":Status.Affected,
                                     "-":Status.Healthy}}
    
class Pedigree:
    
    def __init__(self) -> None:
        self.__pedigree = []

    def read_from_file(self, fname: str) -> None:
        self.__pedigree = []
        page = open(fname, "r").readlines()
        for i in page[1:]:
            t = [k.rstrip() for k in i.split(",")]
            id = int(t[0])
            gender = Gender.Male
            if t[1] == "F":
                gender = Gender.Female
            status = Status.Healthy
            if t[2] == "A":
                status = Status.Affected
            father = None
            mother = None
            if t[3] != "N":
                father = self.__pedigree[int(t[3])]
            if t[4] != 'N':
                mother = self.__pedigree[int(t[4])]
            person = Person(id, gender, status, father, mother)
            assert person not in self.__pedigree, "The same person has been previously added to this pedigree..."
            self.__pedigree.append(person)
            
    def add(self, person: Person) -> None:
        assert person not in self.__pedigree, "The same person has been previously added to this pedigree..."
        self.__pedigree.append(person)
    
    def _pedigree_exists(self, new_pedigree, result_pedigree):
        for existing_pedigree in result_pedigree:
            if all(new_person == existing_person
                   for new_person, existing_person in zip(new_pedigree, existing_pedigree)):
                return True
        return False
    
    def fit(self, model:Models, max_iter:int=100000):
        corrects = {"correct":0,
                    "pedigree":[]}
        if (model["code"] == "AD") or (model['code'] == "AR"):
            for k in range(max_iter):
                for person in self.__pedigree:
                    allele = random.choice(model[person.status])
                    person.somatic.allele = allele
                positive = 0
                for person in self.__pedigree:
                    person.inherit()
                    if model['reverse'][person.somatic.allele] == person.status:
                        positive = positive + 1
                if positive == len(self.__pedigree):
                    corrects['correct'] = corrects['correct'] + 1
                    p = copy.deepcopy(self.__pedigree)
                    if not self._pedigree_exists(p, corrects['pedigree']):
                        corrects['pedigree'].append(p)
        elif model["code"] == "YL":
            pass
        return corrects['correct'] / max_iter, corrects['pedigree']
        
    def generate_random_pedigree(self) -> List[Person]:
        new_pedigree = copy.deepcopy(self.__pedigree)
        for person in new_pedigree:
            person._Person__status = random.choice(list(Status))
        return new_pedigree

    def calculate_fit(self, pedigree: List[Person], model: Models, max_iter: int) -> float:
        correct_count = 0
        for _ in range(max_iter):
            for person in pedigree:
                allele = random.choice(model[person.status])
                person.somatic.allele = allele
            all_correct = True
            for person in pedigree:
                person.inherit()
                if model['reverse'][person.somatic.allele] != person.status:
                    all_correct = False
                    break
            if all_correct:
                correct_count += 1
        return correct_count / max_iter

    def process_chunk(self, chunk_size: int, observed_fit: float, model: Models, max_iter: int) -> int:
        count = 0
        for _ in range(chunk_size):
            random_pedigree = self.generate_random_pedigree()
            random_fit = self.calculate_fit(random_pedigree, model, max_iter)
            if random_fit >= observed_fit:
                count += 1
        return count

    def calc_p_value(self, model: Models, num_simulations: int = 10000, max_iter: int = 1000) -> float:
        observed_fit = self.calculate_fit(self.__pedigree, model, max_iter)
        
        num_cores = multiprocessing.cpu_count()
        chunk_size = num_simulations // num_cores
        
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            futures = [executor.submit(self.process_chunk, chunk_size, observed_fit, model, max_iter) for _ in range(num_cores)]
            count_higher_or_equal = sum(future.result() for future in as_completed(futures))
        
        p_value = (count_higher_or_equal + 1) / (num_simulations + 1)
        return p_value
        
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.__n >= len(self):
            raise StopIteration
        self.__n += 1
        return self.__pedigree[self.__n - 1]
    
    def __getitem__(self, item) -> Person:
        return self.__pedigree[item]
    
    def __len__(self) -> int:
        return len(self.__pedigree)
    
    def __repr__(self) -> str:
        result = ""
        for z in self.__pedigree:
            result = result + str(z) + '\n'
        return result[:-1]
    
if __name__ == "__main__":
    ped = Pedigree()
    ped.read_from_file("ped2.ped")
    result = ped.fit(Models.AutosomalRecessive)
    print(result)
    
    p_value = ped.calc_p_value(Models.AutosomalRecessive)  
    print(p_value)
    
    