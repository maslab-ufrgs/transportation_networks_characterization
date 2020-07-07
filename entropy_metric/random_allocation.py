# -*- coding: utf-8 -*-

import random

def memoize(f):
    memo = {}
    def helper(n,k,i):
        t = (n,k,i) 
        if t not in memo:            
            memo[t] = f(n,k,i)
        return memo[t]
    return helper


def probability_allocation(n, k, i):
    '''
    
    
    Calcula a probabilidade de se alocar aleatoriamente i trips
    (de um conjunto de n trips) na primeira das k rotas que faltam
    ser avaliadas para alocação.

    Essa função é auxiliar da função random_allocation abaixo.
    '''
    if k == 2:
        return 1 / float(n+1)
        # Observe que se há duas rotas que restam ser consideradas,
        # então há apenas um *grau de liberdade* nessa escolha, uma
        # vez que a alocação deve satisfazer a restrição da demanda,
        # isto é, a soma das alocações é precisamente a demanda.
    else:
        Probability = 1 - n / float(n + k - 1)
        for j in range(1,k-1):
            Probability *= (1- i/float(n+j))
        return Probability

# Dando a propriedade memoize para a função probability allocation
probability_allocation = memoize( probability_allocation )
    
def random_allocation(n,k,g=1):
    '''
    Sorteia uma tupla de *k* entradas inteiras e maiores ou iguais
    a zero cujas entradas totalizam *n*.
    '''
    if k == 1:
        # Esse é o caso base do processo de recursão.
        # Quando há apenas uma rota que falta ser considerada,
        # então basta retornar o número de trips que não foi
        # alocado em nenhuma das demais rotas.
        return [n,]
    elif n == 0:
        # Esse é o caso do encerramento prematuro da aleatoriedade,
        # pois todas as trips já foram alocadas nas rotas previamente
        # consideradas. resta então retornar uma lista de k entradas iguais a 0
        return k * [0,]
    else:
        # Aqui ocorre a chamada recursiva da função
        q = random.random()
        acc_probability = 0
        for i in range(0,n+1):
            acc_probability += probability_allocation(n,k,i)
            if q < acc_probability:
                # Observe que é preciso subtrair a demanda
                # na chamada recursiva, assim como subtrair 1
                # do número de rotas que será considerado na
                # próxima chamada recursiva.
                return [i,] + random_allocation(n-i,k-1)
        return [n,] + random_allocation(0,k-1)

def random_allocation_extended(n,k,g=1):
    
    '''  
    Sorteia  uma tupla  de *k*  entradas
    inteiras e maiores  ou iguais a zero
    cujas entradas totalizam *n*.
    '''
    
    r_allocation = random_allocation(n,k)
    # print r_allocation

    extended_allocation = list()

    i = 0
    for a in r_allocation:
        extended_allocation += a * [i,]
        i += 1
        
    return extended_allocation
    
    
if __name__ == '__main__':
    print("Random_allocation_script")