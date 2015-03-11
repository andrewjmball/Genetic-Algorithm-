
from __future__ import division
import numpy
import cmath
from pyevolve import G1DList, GSimpleGA, Selectors
from pyevolve import Initializators, Mutators, Consts
from pyevolve import DBAdapters 

from pyevolve import Scaling

A = 0.067 #EffectiveMassinGaAs
B = 0.57  #EffectiveMassinDiamond
C = 1    #EffectiveMassinVacuum

def eval_func1(ind):
    score = 0.0
    if ind[0] > 5 and ind[1] > ind[0] and ind[2] > 1+ ind[1]:
            score = -1  
    return score

def eval_func2(ind):
    data_position = [0,0,0,0]
    data_position.pop()
    data_position.pop()
    data_position.pop()
    data_position.append(ind[0])
    data_position.append(ind[1])
    data_position.append(ind[2])
    
    data_V = [0]
    data_V.append(4.046)   #setting potential voltages
    data_V.append(4.046)
    data_V.append(0)
    data_V.append(0)
    data_V.append(4)
    data_V.append(4)
    data_V.append(-4)
    
    
    input_ninterfaces = 4                                                                 
    data_EffMass = [A,B,B,A,A,C,C,A]           # initial effective masses before vacuum, each interface has 2 values, for left and right of the interface
# intial positions before vacuum, each interface has 1 value

                                                          
    def calculatert(E):                                       # correctly working code below tested on milestone
        matrix_list = []
        for i in range(int(input_ninterfaces)):                              
            Mminus = data_EffMass[2*i]
            Mplus = data_EffMass[2*i+1]
            Kminus = cmath.sqrt((E-data_V[2*i])*Mminus/13.6057)
            Kplus = cmath.sqrt((E-data_V[2*i+1])*Mplus/13.6057)
            
            #below creates the matricies for each interface  
            A_array = numpy.array([[numpy.exp(1j*Kplus*data_position[i]),numpy.exp(-1j*Kplus*data_position[i])],[(Kplus/Mplus)*numpy.exp(1j*Kplus*data_position[i]),-Kplus/Mplus*numpy.exp(-1j*Kplus*data_position[i])]])
            A = numpy.asmatrix(A_array)
            B = A.I
            C_array = numpy.array([[numpy.exp(1j*Kminus*data_position[i]),numpy.exp(-1j*Kminus*data_position[i])],[(Kminus/Mminus)*numpy.exp(1j*Kminus*data_position[i]),-Kminus/Mminus*numpy.exp(-1j*Kminus*data_position[i])]])
            C = numpy.asmatrix(C_array)
            D = B*C
            matrix_list.append(D)  
    
        end_matrix = numpy.identity(2) 
        for i in range(int(input_ninterfaces)): 
            end_matrix = matrix_list[i] * end_matrix                         #this is the matrix multiplication for all interfaces
            
        Littler = -end_matrix[1,0]/end_matrix[1,1]                           #this code creates the T and R coefficients 
        Bigr = (abs(Littler))**2
        Littlet = end_matrix[0,0] - ((end_matrix[0,1]*end_matrix[1,0])/ end_matrix[1,1])
        Bigt = (abs(Littlet))**2
        return Bigr, Bigt
   
    score = calculatert(3.730)[0] #minimize reflection at 3.730 (maximize transmission)
    return score  

def eval_func3(ind):
    data_position = [0,0,0,0]
    data_position.pop()
    data_position.pop()
    data_position.pop()
    data_position.append(ind[0])
    data_position.append(ind[1])
    data_position.append(ind[2])
    
    data_V = [0]
    data_V.append(4.046)
    data_V.append(4.046)
    data_V.append(0)
    data_V.append(0)
    data_V.append(4)
    data_V.append(4)
    data_V.append(-4)
    
    
    input_ninterfaces = 4                                                                 
    data_EffMass = [A,B,B,A,A,C,C,A]           # initial effective masses before vacuum, each interface has 2 values, for left and right of the interface
# intial positions before vacuum, each interface has 1 value

                                                          
    def calculatert(E):                                       # correctly working code below tested on milestone
        matrix_list = []
        for i in range(int(input_ninterfaces)):                              
            Mminus = data_EffMass[2*i]
            Mplus = data_EffMass[2*i+1]
            Kminus = cmath.sqrt((E-data_V[2*i])*Mminus/13.6057)
            Kplus = cmath.sqrt((E-data_V[2*i+1])*Mplus/13.6057)
            
            #below creates the matricies for each interface  
            A_array = numpy.array([[numpy.exp(1j*Kplus*data_position[i]),numpy.exp(-1j*Kplus*data_position[i])],[(Kplus/Mplus)*numpy.exp(1j*Kplus*data_position[i]),-Kplus/Mplus*numpy.exp(-1j*Kplus*data_position[i])]])
            A = numpy.asmatrix(A_array)
            B = A.I
            C_array = numpy.array([[numpy.exp(1j*Kminus*data_position[i]),numpy.exp(-1j*Kminus*data_position[i])],[(Kminus/Mminus)*numpy.exp(1j*Kminus*data_position[i]),-Kminus/Mminus*numpy.exp(-1j*Kminus*data_position[i])]])
            C = numpy.asmatrix(C_array)
            D = B*C
            matrix_list.append(D)  
    
        end_matrix = numpy.identity(2) 
        for i in range(int(input_ninterfaces)): 
            end_matrix = matrix_list[i] * end_matrix                         #this is the matrix multiplication for all interfaces
            
        Littler = -end_matrix[1,0]/end_matrix[1,1]                           #this code creates the T and R coefficients 
        Bigr = (abs(Littler))**2
        Littlet = end_matrix[0,0] - ((end_matrix[0,1]*end_matrix[1,0])/ end_matrix[1,1])
        Bigt = (abs(Littlet))**2
        return Bigr, Bigt
   
    score = calculatert(2.986)[1] #minimize tranmission at 2.986 (maximize reflection)
    return score  
    

               
# Genome instance
genome = G1DList.G1DList(3)
genome.setParams(rangemin=0, rangemax=500, bestRawScore=0.00, roundDecimal=2)
genome.initializator.set(Initializators.G1DListInitializatorReal)
genome.mutator.set(Mutators.G1DListMutatorRealGaussian)

# The evaluator function (objective function)
genome.evaluator.set(eval_func2)
genome.evaluator.add(eval_func1)
genome.evaluator.add(eval_func3)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.selector.set(Selectors.GRouletteWheel)

#logging and committing to DB
sqlite_adapter = DBAdapters.DBSQLite(identify="ex1")
ga.setDBAdapter(sqlite_adapter)

pop = ga.getPopulation()
pop.scaleMethod.set(Scaling.SigmaTruncScaling)

ga.minimax = Consts.minimaxType["minimize"]
ga.setGenerations(500)
ga.setMutationRate(0.05)

#hold in database 
sqlite_adapter = DBAdapters.DBSQLite(identify="ex1")
ga.setDBAdapter(sqlite_adapter)

# Do the evolution, with stats dump
# frequency of 10 generations
ga.evolve(freq_stats=100)

# Best individual
best = ga.bestIndividual()
print "\nBest individual score: %.2f" % (best.score,)
print best