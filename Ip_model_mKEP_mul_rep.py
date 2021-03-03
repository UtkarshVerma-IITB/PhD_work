# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 03:17:47 2020

@author: Dr. Utkarsh Verma
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 01:03:40 2020

@author: Dr. Utkarsh Verma
"""

from pulp import *
import random
import numpy as np
import time
import matplotlib.pyplot as plt

rng1 = np.random.RandomState(11)


# =============================================================================
# n=[1,2,3,4,5]
# nc1=[1,2,3]
# nc2=[4,5]
# arcs={(1,2),(1,4),(2,1),(2,3),(3,2),(3,5),(4,1),(5,3)}
# c={(1,2): 1,
#    (1,4): 2,
#    (2,1): 1,
#    (2,3): 1,
#    (3,2): 1,
#    (3,5): 2,
#    (4,1): 2,
#    (5,3): 2}
# M=[1,2,3]
# =============================================================================


def IP_model_int(n,nc1,nc2,arcs,c,M,b1,b2):
        
    #print ('arcs',arcs,'c',c)
    prob1= LpProblem("Int_model_for_KEP", LpMaximize)
    x = LpVariable.dicts('Route var',(arcs,M),0,1, cat='Integer')
    y1= LpVariable('var1',0, cat='Continous')
    y2= LpVariable('var2',0, cat ='Continous')
    
# =============================================================================
#     for i in arcs:
#         if i==(1,2):
#             for j in M:
#                 print ('x[i]',x[i][j],i[0],i[1])
#     
# =============================================================================
    # =============================================================================
    # for i in arcs:
    #     print ('i',i,i[0],i[1])
    #     print ('c',c[i])
    #     for j in M:
    #         print ('j',j)
    #         print ('x',x[i][j])
    # =============================================================================
    
    prob1 += lpSum([x[i][j]*c[i] for i in arcs for j in M ])#-100*(-lpSum([x[i][l] for i in arcs for l in M for k in nc1 if i[0]==k]) + b1) -100*(-lpSum([x[i][l] for i in arcs for l in M for k in nc2 if i[0]==k]) + b2) 
    
# =============================================================================
#     for j in nc1:
#         for k in arcs:
#             if k[0]==j:
#                 print ('k',k,'k[0]',k[0])
#     
# =============================================================================
    for k in n:
        for j in M:
            prob1 += lpSum([x[i][j] for i in arcs if i[0]==k]) - lpSum([x[i][j] for i in arcs if i[1]==k]) ==0
    
    for k in n:
        prob1 += lpSum([x[i][j] for i in arcs for j in M if i[0]==k]) <=1
        
    for l in M:
        prob1 += lpSum([x[i][l] for i in arcs]) <=3
        
    for l in M:
        prob1 += lpSum(x[i][l] for i in arcs if i[0] in nc1 if i[1] in nc1) <=3
    
    prob1 += lpSum([x[i][l] for i in arcs for l in M for k in nc1 if i[0]==k])   >=b1
    prob1 += lpSum([x[i][l] for i in arcs for l in M for k in nc2 if i[0]==k])   >=b2
        
    prob1.writeLP('Matching.LP')
    
    prob1.solve()
    #solver = pl.CPLEX_PY()
    #prob1.solve(solver)

    print ('Status',LpStatus[prob1.status])
    print ('Objective',value(prob1.objective))
    #print ('x',prob1.variables())

    z_name=[]
    z_value=[]
    for i in prob1.variables():
        if i.value()!=0:
            z_name.append(i.name)
            z_value.append(i.value())
            print (i,i.value())
    
    return z_name,z_value,value(prob1.objective)

#Ex,GS=IP_model_int(n,nc1,nc2,arcs,c,M)

#print ('Ex',Ex,'Global Solution',GS)

def IP_model_ind(n,arcs,c,M):
        
    #print ('arcs',arcs,'c',c)
    prob1= LpProblem("Int_model_for_KEP", LpMaximize)
    x = LpVariable.dicts('Route var',(arcs,M),0,1, cat='Integer')
    
# =============================================================================
#     for i in arcs:
#         if i==(1,2):
#             for j in M:
#                 print ('x[i]',x[i][j],i[0],i[1])
#     
# =============================================================================
    # =============================================================================
    # for i in arcs:
    #     print ('i',i,i[0],i[1])
    #     print ('c',c[i])
    #     for j in M:
    #         print ('j',j)
    #         print ('x',x[i][j])
    # =============================================================================
    
    prob1 += lpSum([x[i][j]*c[i] for i in arcs for j in M])
    
# =============================================================================
#     for j in nc1:
#         for k in arcs:
#             if k[0]==j:
#                 print ('k',k,'k[0]',k[0])
#     
# =============================================================================
    for k in n:
        for j in M:
            prob1 += lpSum([x[i][j] for i in arcs if i[0]==k]) - lpSum([x[i][j] for i in arcs if i[1]==k]) ==0
    
    for k in n:
        prob1 += lpSum([x[i][j] for i in arcs for j in M if i[0]==k]) <=1
        
    for l in M:
        prob1 += lpSum([x[i][l] for i in arcs]) <=3
                
    prob1.writeLP('Matching.LP')
    
    prob1.solve()
    #solver = pl.CPLEX_PY()
    #prob1.solve(solver)
    
    print ('Status GS',LpStatus[prob1.status])
    print ('Objective GS',value(prob1.objective))
    #print ('x',prob1.variables())
    
    z_name=[]
    z_value=[]
    for i in prob1.variables():
        if i.value()!=0:
            z_name.append(i.name)
            z_value.append(i.value())
            print (i,i.value())
    
    return z_name,z_value,value(prob1.objective)

def generating_bg_ASTRA(n):
    x=[4.0/211,44.0/211,59.0/211,60.0/211,102.0/211,104.0/211,115.0/211,119.0/211,120.0/211,121.0/211,122.0/211,123.0/211,
       156.0/211,200.0/211,210.0/211,1]
    generated_group=[]
    generated_rec=[]
    generated_don=[]
    for i in range(n):
        random_num=rng1.uniform(0,1)
        if random_num>=0 and random_num <=x[0]:
            #print ('generated group is A,A')
            generated_group.append('A,A')
            generated_rec.append('A')
            generated_don.append('A')
        elif random_num>=x[0] and random_num <=x[1]:
            
            #print ('generated group is A,B')
            generated_group.append('A,B')
            generated_rec.append('A')
            generated_don.append('B')
            
        elif random_num>=x[1] and random_num <=x[2]:
            
            #print 'generated group is A,AB'
            generated_group.append('A,AB')
            generated_rec.append('A')
            generated_don.append('AB')
            
        elif random_num>=x[2] and random_num <=x[3]:
            
            #print 'generated group is A,O'
            generated_group.append('A,O')
            generated_rec.append('A')
            generated_don.append('O')
            
        elif random_num>=x[3] and random_num <=x[4]:
            
            #print 'generated group is B,A'
            generated_group.append('B,A')
            generated_rec.append('B')
            generated_don.append('A')
            
        elif random_num>=x[4] and random_num <=x[5]:
            
            #print 'generated group is B,B'
            generated_group.append('B,B')
            generated_rec.append('B')
            generated_don.append('B')
            
        elif random_num>=x[5] and random_num <=x[6]:
            #print 'generated group is B,AB'
            generated_group.append('B,AB')
            generated_rec.append('B')
            generated_don.append('A')
            
        elif random_num>=x[6] and random_num <=x[7]:
            
            #print 'generated group is B,O'
            generated_group.append('B,O')
            generated_rec.append('B')
            generated_don.append('O')
            
        elif random_num>=x[7] and random_num <=x[8]:
            
            #print 'generated group is AB,A'
            generated_group.append('AB,A')
            generated_rec.append('AB')
            generated_don.append('A')
            
        elif random_num>=x[8] and random_num <=x[9]:
            
            #print 'generated group is AB,B'
            generated_group.append('AB,B')
            generated_rec.append('AB')
            generated_don.append('B')
            
        elif random_num>=x[9] and random_num <=x[10]:
            
            #print 'generated group is AB,AB'
            generated_group.append('AB,AB')
            generated_rec.append('AB')
            generated_don.append('AB')
            
        elif random_num>=x[10] and random_num <=x[11]:
            
            #print 'generated group is A,A'
            generated_group.append('AB,O')
            generated_rec.append('AB')
            generated_don.append('O')
            
        elif random_num>=x[11] and random_num <=x[12]:
            
            #print 'generated group is O,A'
            generated_group.append('O,A')
            generated_rec.append('O')
            generated_don.append('A')
            
        elif random_num>=x[12] and random_num <=x[13]:
            
            #print 'generated group is O,B'
            generated_group.append('O,B')
            generated_rec.append('O')
            generated_don.append('B')
            
        elif random_num>=x[13] and random_num <=x[14]:
            
            #print 'generated group is O,AB'
            generated_group.append('O,AB')
            generated_rec.append('O')
            generated_don.append('AB')
        elif random_num>=x[14] and random_num <=1:
            
            #print 'generated group is O,O'
            generated_group.append('O,O')
            generated_rec.append('O')
            generated_don.append('O')
        else:
            print ('No pair is generated')
    return generated_rec,generated_don



def generating_bg_saidman(n):
    x=[0.23174596, 0.39412218, 0.4628661,  0.4814, 0.64377622, 0.75754751, 
       0.80571395, 0.8187, 0.88744392, 0.93561036, 0.9560022,  0.9615,
       0.9800339,  0.99301995, 0.99851775, 1]
# =============================================================================
#     x=[4.0/211,44.0/211,59.0/211,60.0/211,102.0/211,104.0/211,115.0/211,119.0/211,120.0/211,121.0/211,122.0/211,123.0/211,
#        156.0/211,200.0/211,210.0/211,1]
# =============================================================================
    y=[0.37,0.23,0.32,0.08] #O,A,B,AB
    y1=[]
    for i in y:
        for j in y:
            y1.append(i*j)
            
    x=np.cumsum(y1)
    generated_group=[]
    generated_rec=[]
    generated_don=[]
    for i in range(n):
        random_num=rng1.uniform(0,1)
        if random_num>=0 and random_num <=x[0]:
            #print ('generated group is A,A')
            generated_group.append('O,O')
            generated_rec.append('O')
            generated_don.append('O')
        elif random_num>=x[0] and random_num <=x[1]:
            
            #print ('generated group is A,B')
            generated_group.append('O,A')
            generated_rec.append('O')
            generated_don.append('A')
            
        elif random_num>=x[1] and random_num <=x[2]:
            
            #print 'generated group is A,AB'
            generated_group.append('O,B')
            generated_rec.append('O')
            generated_don.append('B')
            
        elif random_num>=x[2] and random_num <=x[3]:
            
            #print 'generated group is A,O'
            generated_group.append('O,AB')
            generated_rec.append('O')
            generated_don.append('AB')
            
        elif random_num>=x[3] and random_num <=x[4]:
            
            #print 'generated group is B,A'
            generated_group.append('A,O')
            generated_rec.append('A')
            generated_don.append('O')
            
        elif random_num>=x[4] and random_num <=x[5]:
            
            #print 'generated group is B,B'
            generated_group.append('A,A')
            generated_rec.append('A')
            generated_don.append('A')
            
        elif random_num>=x[5] and random_num <=x[6]:
            #print 'generated group is B,AB'
            generated_group.append('A,B')
            generated_rec.append('A')
            generated_don.append('B')
            
        elif random_num>=x[6] and random_num <=x[7]:
            
            #print 'generated group is B,O'
            generated_group.append('A,AB')
            generated_rec.append('A')
            generated_don.append('AB')
            
        elif random_num>=x[7] and random_num <=x[8]:
            
            #print 'generated group is AB,A'
            generated_group.append('B,O')
            generated_rec.append('B')
            generated_don.append('O')
            
        elif random_num>=x[8] and random_num <=x[9]:
            
            #print 'generated group is AB,B'
            generated_group.append('B,A')
            generated_rec.append('B')
            generated_don.append('A')
            
        elif random_num>=x[9] and random_num <=x[10]:
            
            #print 'generated group is AB,AB'
            generated_group.append('B,B')
            generated_rec.append('B')
            generated_don.append('B')
            
        elif random_num>=x[10] and random_num <=x[11]:
            
            #print 'generated group is A,A'
            generated_group.append('B,AB')
            generated_rec.append('B')
            generated_don.append('AB')
            
        elif random_num>=x[11] and random_num <=x[12]:
            
            #print 'generated group is O,A'
            generated_group.append('AB,O')
            generated_rec.append('AB')
            generated_don.append('O')
            
        elif random_num>=x[12] and random_num <=x[13]:
            
            #print 'generated group is O,B'
            generated_group.append('AB,A')
            generated_rec.append('AB')
            generated_don.append('A')
            
        elif random_num>=x[13] and random_num <=x[14]:
            
            #print 'generated group is O,AB'
            generated_group.append('AB,B')
            generated_rec.append('AB')
            generated_don.append('B')
        elif random_num>=x[14] and random_num <=1:
            
            #print 'generated group is O,O'
            generated_group.append('AB,AB')
            generated_rec.append('AB')
            generated_don.append('AB')
        else:
            print ('No pair is generated')
    return generated_rec,generated_don

def compatible_edges(recipient_bg, donor_bg):
    Edges=[]
    for i in range(0,len(donor_bg)):
        for j in range(0,len(recipient_bg)):
            if ((donor_bg[i]==recipient_bg[j] or donor_bg[i]=='O' or recipient_bg[j]=='AB') and (donor_bg[i]!="Matched" and recipient_bg[i]!="Matched" and donor_bg[j]!="Matched" and recipient_bg[j]!="Matched")):
                Edges.append((i,j))
                continue

    return Edges
  
#Function for creating age difference
def generate_rec_age():
    r=rng1.uniform(0,1)
    if r <= 4/284: #Freq = 4
        r1=rng1.randint(10,20)
    elif r >=4/284 and r <= 51/284: #Frequency = 47
        r1=rng1.randint(20,30)
    elif r >=51/284 and r <= 121/284: #Freq = 70
        r1=rng1.randint(30,40)
    elif r >=121/284 and r <= 194/284: #freq = 73 
        r1=rng1.randint(40,50)
    elif r >=194/284 and r <= 263/284: #freq = 69
        r1=rng1.randint(50,60)
    elif r >=263/284: #freq = 21
        r1=rng1.randint(60,70)
    return r1

def generate_don_age():
    r=rng1.uniform(0,1)
    if r <= 15/283: #Freq = 15
        r1=rng1.randint(20,30)
    elif r >=15/283 and r <= 54/283: #Frequency = 39
        r1=rng1.randint(30,40)
    elif r >=54/283 and r <= 163/283: #Freq = 109
        r1=rng1.randint(40,50)
    elif r >=163/283 and r <= 251/283: #freq = 88 
        r1=rng1.randint(50,60)
    elif r >=251/283 and r<=1: #freq = 32
        r1=rng1.randint(60,70)
    return r1
    
def HLA6_type():
    HLA6=[]
    
    a1=[17/80,32/80,42/80,58/80,72/80,75/80,77/80,79/80,1] #HLA alleles [1,2,3,11,24,26,30,31,33]
    a2=[3/72,6/72,7/72,16/72,24/72,25/72,35/72,42/72,64/72,1] #HLA alleles [2,3,7,11,24,26,29,32,33,68]
    b1=[26/80,32/80,35/80,43//80,44//80,52/80,54/80,56/80,63/80,72/80,73/80,75/80,1] #HLA alleles [7,8,13,15,33,35,37,39,40,44,50,51,52]
    b2=[2/78,4/78,6/78,9/78,16/78,20/78,33/78,44/78,46/78,60/78,66/78,67/78,72/78,77/78,1] #HLA alleles [7,8,15,27,35,37,40,44,51,52,55,56,57,58,60]
    dr1=[6/81,15/81,28/81,29/81,32/81,50/81,51/81,58/81,60/81,62/81,66/81,71/81,1] #HLA alleles [1,3,4,5,6,7,9,10,11,12,13,14,15]
    dr2=[2/78,6/78,18/78,19/78,28/78,33/78,35/78,45/78,52/78,1] #HLA alleles [3,4,7,9,10,11,12,13,14,15]
    #Generation of HLA A1 alleles
    random_a1=rng1.uniform(0,1)
    if random_a1<=a1[0]:
        HLA6.append("1")
    elif random_a1 >=a1[0] and random_a1<=a1[1]: 
        HLA6.append("2")
    elif random_a1 >=a1[1] and random_a1<=a1[2]: 
        HLA6.append("3")
    elif random_a1 >=a1[2] and random_a1<=a1[3]: 
        HLA6.append("11")
    elif random_a1 >=a1[3] and random_a1<=a1[4]: 
        HLA6.append("24")
    elif random_a1 >=a1[4] and random_a1<=a1[5]: 
        HLA6.append("26")
    elif random_a1 >=a1[5] and random_a1<=a1[6]: 
        HLA6.append("30")
    elif random_a1 >=a1[6] and random_a1<=a1[7]: 
        HLA6.append("31")
    elif random_a1 >=a1[7] and random_a1<=a1[8]: 
        HLA6.append("33")
    #Generation of HLA A2 alleles
    random_a2=rng1.uniform(0,1)
    if  random_a2<=a2[0]:
        HLA6.append("2")
    elif random_a2 >=a2[0] and random_a2<=a2[1]: 
        HLA6.append("3")
    elif random_a2 >=a2[1] and random_a2<=a1[2]: 
        HLA6.append("7")
    elif random_a2 >=a2[2] and random_a2<=a2[3]: 
        HLA6.append("11")
    elif random_a2 >=a2[3] and random_a2<=a2[4]: 
        HLA6.append("24")
    elif random_a2 >=a2[4] and random_a2<=a2[5]: 
        HLA6.append("26")
    elif random_a2 >=a2[5] and random_a2<=a2[6]: 
        HLA6.append("29")
    elif random_a2 >=a2[6] and random_a2<=a2[7]: 
        HLA6.append("32")
    elif random_a2 >=a2[7] and random_a2<=a2[8]: 
        HLA6.append("33")
    elif random_a2 >=a2[8] and random_a2<=a2[9]: 
        HLA6.append("68")
    #Generation of HLA B1 alleles    
    random_b1=rng1.uniform(0,1)
    if random_b1 >=0 and random_b1<=b1[0]:
        HLA6.append("7")
    elif random_b1 >=b1[0] and random_b1<=b1[1]: 
        HLA6.append("8")
    elif random_b1 >=b1[1] and random_b1<=b1[2]: 
        HLA6.append("13")
    elif random_b1 >=b1[2] and random_b1<=b1[3]: 
        HLA6.append("15")
    elif random_b1 >=b1[3] and random_b1<=b1[4]: 
        HLA6.append("33")
    elif random_b1 >=b1[4] and random_b1<=b1[5]: 
        HLA6.append("35")
    elif random_b1 >=b1[5] and random_b1<=b1[6]: 
        HLA6.append("37")
    elif random_b1 >=b1[6] and random_b1<=b1[7]: 
        HLA6.append("39")
    elif random_b1 >=b1[7] and random_b1<=b1[8]: 
        HLA6.append("40")
    elif random_b1 >=b1[8] and random_b1<=b1[9]: 
        HLA6.append("44")
    elif random_b1 >=b1[9] and random_b1<=b1[10]: 
        HLA6.append("50")
    elif random_b1 >=b1[10] and random_b1<=b1[11]: 
        HLA6.append("51")
    elif random_b1 >=b1[11] and random_b1<=b1[12]: 
        HLA6.append("52")
    #Generation of HLA B2 alleles    
    random_b2=rng1.uniform(0,1)
    if random_b2 >=0 and random_b2<=b1[0]:
        HLA6.append("7")
    elif random_b2 >=b2[0] and random_b2<=b2[1]: 
        HLA6.append("8")
    elif random_b2 >=b2[1] and random_b2<=b2[2]: 
        HLA6.append("15")
    elif random_b2 >=b2[2] and random_b2<=b2[3]: 
        HLA6.append("27")
    elif random_b2 >=b2[3] and random_b2<=b2[4]: 
        HLA6.append("35")
    elif random_b2 >=b2[4] and random_b2<=b2[5]: 
        HLA6.append("37")
    elif random_b2 >=b2[5] and random_b2<=b2[6]: 
        HLA6.append("40")
    elif random_b2 >=b2[6] and random_b2<=b2[7]: 
        HLA6.append("44")
    elif random_b2 >=b2[7] and random_b2<=b2[8]: 
        HLA6.append("51")
    elif random_b2 >=b2[8] and random_b2<=b2[9]: 
        HLA6.append("52")
    elif random_b2 >=b2[9] and random_b2<=b2[10]: 
        HLA6.append("55")
    elif random_b2 >=b2[10] and random_b2<=b2[11]: 
        HLA6.append("56")
    elif random_b2 >=b2[11] and random_b2<=b2[12]: 
        HLA6.append("57")
    elif random_b2 >=b2[12] and random_b2<=b2[13]: 
        HLA6.append("58")
    elif random_b2 >=b2[13] and random_b2<=b2[14]: 
        HLA6.append("60")
    #Generation of HLA DR1 alleles    
    random_dr1=rng1.uniform(0,1)
    if random_dr1 >=0 and random_dr1<=dr1[0]:
        HLA6.append("1")
    elif random_dr1 >=dr1[0] and random_dr1<=dr1[1]: 
        HLA6.append("3")
    elif random_dr1 >=dr1[1] and random_dr1<=dr1[2]: 
        HLA6.append("4")
    elif random_dr1 >=dr1[2] and random_dr1<=dr1[3]: 
        HLA6.append("5")
    elif random_dr1 >=dr1[3] and random_dr1<=dr1[4]: 
        HLA6.append("6")
    elif random_dr1 >=dr1[4] and random_dr1<=dr1[5]: 
        HLA6.append("7")
    elif random_dr1 >=dr1[5] and random_dr1<=dr1[6]: 
        HLA6.append("9")
    elif random_dr1 >=dr1[6] and random_dr1<=dr1[7]: 
        HLA6.append("10")
    elif random_dr1 >=dr1[7] and random_dr1<=dr1[8]: 
        HLA6.append("11")
    elif random_dr1 >=dr1[8] and random_dr1<=dr1[9]: 
        HLA6.append("12")
    elif random_dr1 >=dr1[9] and random_dr1<=dr1[10]: 
        HLA6.append("13")
    elif random_dr1 >=dr1[10] and random_dr1<=dr1[11]: 
        HLA6.append("14")
    elif random_dr1 >=dr1[11] and random_dr1<=dr1[12]: 
        HLA6.append("15")
    #Generation of HLA DR2 alleles    
    random_dr2=rng1.uniform(0,1)
    if random_dr2 >=0 and random_dr2<=dr2[0]:
        HLA6.append("3")
    elif random_dr2 >=dr2[0] and random_dr2<=dr2[1]: 
        HLA6.append("4")
    elif random_dr2 >=dr2[1] and random_dr2<=dr2[2]: 
        HLA6.append("7")
    elif random_dr2 >=dr2[2] and random_dr2<=dr2[3]: 
        HLA6.append("9")
    elif random_dr2 >=dr2[3] and random_dr2<=dr2[4]: 
        HLA6.append("10")
    elif random_dr2 >=dr2[4] and random_dr2<=dr2[5]: 
        HLA6.append("11")
    elif random_dr2 >=dr2[5] and random_dr2<=dr2[6]: 
        HLA6.append("12")
    elif random_dr2 >=dr2[6] and random_dr2<=dr2[7]: 
        HLA6.append("13")
    elif random_dr2 >=dr2[7] and random_dr2<=dr2[8]: 
        HLA6.append("14")
    elif random_dr2 >=dr2[8] and random_dr2<=dr2[9]: 
        HLA6.append("15")
        
    return HLA6

def create_score(donor,recipient):
    donor_age=generate_don_age()
    rec_age=generate_don_age()
    donor_HLA=HLA6_type()
    rec_HLA=HLA6_type()
    #print ('donor HLA',donor_HLA)
    #print ('rec HLA',rec_HLA)
    total_score=0
    #Age score calculation
    #age_score=0
    if abs(donor_age - rec_age) <= 40:
        age_score=-(5/4)*abs(donor_age - rec_age) + 50
    else:
        age_score=0
        
# =============================================================================
#     if donor_age - rec_age <= -20 or donor_age - rec_age >= 20 :
#         age_score=15
#     if (donor_age - rec_age >= -20 and donor_age - rec_age <= -10) or (donor_age - rec_age >= 10 and donor_age - rec_age <= 20) :
#         age_score=30
#     if (donor_age - rec_age >= -10 and donor_age - rec_age <= 0) or (donor_age - rec_age >= 0 and donor_age - rec_age <= 10) :
#         age_score=50
# =============================================================================
        
    #HLA score calculation 
    HLA_score=0
    count=0
    for i in range(6):
        if donor_HLA[i]==rec_HLA[i]:
            count=count+1
    if count ==0:
        HLA_score=10
    if count ==1:
        HLA_score=25
    if count ==2:
        HLA_score=40
    if count ==3:
        HLA_score=55
    if count ==4:
        HLA_score=70
    if count ==5:
        HLA_score=85
    if count ==6:
        HLA_score=100
    
    total_score=age_score+HLA_score
    
    return total_score

print ('Total score',create_score('D1','R1'))

def count_bg(BG_list,BG):
    count=0
    for i in BG_list:
        if i==BG:
            count=count+1
    return count

def managing_pairs(arcs,cost): #To remove compatible edges which has score less than individual score
    CE=[]
    for i in arcs:
        if i[0]==i[1]:
            CE.append(i)
    to_be_removed=[]        
    for j in arcs:
        for i in CE:
            if j[1]==i[1] and cost[j] <= cost[i] and j!=i:
                #print ('edges to be removed',j,cost[j],cost[i])
                to_be_removed.append(j)
                
    for i in to_be_removed:
        arcs.remove(i)
        cost.pop(i)
    
    return arcs, cost

Avg_IS_C1=[]
Avg_IS_C2=[]
Avg_GS_C1=[]
Avg_GS_C2=[]
Avg_IS_C1_score=[]
Avg_IS_C2_score=[]
Avg_GS_C1_score=[]
Avg_GS_C2_score=[]

Avg_IS_C1_matched=[]
Avg_IS_C2_matched=[]
Avg_GS_C1_matched=[]
Avg_GS_C2_matched=[]

Avg_BG_IS_C1=[]
Avg_BG_IS_C2=[]
Avg_BG_GS_C1=[]
Avg_BG_GS_C2=[]

rep=2
runs=2

for b in range(rep):
    print ('Replication',b)
    IS_c1=[]
    IS_c2=[]
    Num_matched_c1=[]
    Num_matched_c2=[]
    GS=[]
    GS_c1=[]
    GS_c2=[]
    GS_c1_score=[]
    GS_c2_score=[]
    c1_rec_bg=[]
    c1_don_bg=[]
    c2_rec_bg=[]
    c2_don_bg=[]
    n_rec=[]
    n_don=[]
    p1=0.2 #dropout probability for country 1
    p2=0.2 #dropout probability for country 2
    
    GR_BG_Rec_country1=[]
    GR_BG_Rec_country2=[]
    
    GR_BG_Don_country1=[]
    GR_BG_Don_country2=[]
    
    BG_rec_c1=[]
    BG_rec_c2=[]
    
    GS_BG_rec_c1=[]
    GS_BG_rec_c2=[]
    
    #start_time = time.clock()
    #z=3
    print ('New Program')
    for a in range(runs):
        print ('Round',a)
        a1=rng1.randint(5,10)#Random number of pairs generated for C1
        print ('number of pairs in c1',a1)
        a2=rng1.randint(5,10)#Random number of pairs generated for C2
        print ('number of pairs in c2',a2)
        nc1_rec, nc1_don=generating_bg_ASTRA(a1) #BG generation of those pairs
        nc2_rec, nc2_don=generating_bg_ASTRA(a2)
        
        print ('Country 1 pairs DR',nc1_don,nc1_rec)
        print ('Country 2 pairs DR',nc2_don,nc2_rec)
        
        c1_rec_bg=c1_rec_bg+nc1_rec #Number of pairs considered for matching
        c1_don_bg=c1_don_bg+nc1_don
        c2_rec_bg=c2_rec_bg+nc2_rec
        c2_don_bg=c2_don_bg+nc2_don
        
        nc1=[i for i in range(len(c1_rec_bg))] #Number of nodes in C1
        arcs_c1=compatible_edges(c1_rec_bg,c1_don_bg)
        cost_c1={(i,j):create_score(i,j) for (i,j) in arcs_c1}
        M_c1=[i for i in range(1,int((len(nc1)+1)/2)+1)]
        arcs_c1_new, cost_c1_new = managing_pairs(arcs_c1,cost_c1)
        Ex1_c1_name,Ex1_c1_value, IS1 =IP_model_ind(nc1,arcs_c1_new,cost_c1_new,M_c1)
        print ('IS1 before',IS1)
        
        if str(IS1)=='None':
            IS1=0
        print ('IS1 after',IS1)
        IS_c1.append(IS1)
        
        Node_c1=[]
        for i in range(len(Ex1_c1_name)):
            if Ex1_c1_name[i]!='__dummy':
                print (Ex1_c1_name[i].split('_')[2].split('(')[1].split(',')[0],"=",Ex1_c1_value[i])
                Node_c1.append(int(Ex1_c1_name[i].split('_')[2].split('(')[1].split(',')[0]))
                BG_rec_c1.append(c1_rec_bg[int(Ex1_c1_name[i].split('_')[2].split('(')[1].split(',')[0])])
        Num_matched_c1.append(len(Node_c1))
    
        nc2=[i for i in range(len(c2_rec_bg))] #Number of nodes in C2
        arcs_c2=compatible_edges(c2_rec_bg,c2_don_bg)
        cost_c2={(i,j):create_score(i,j) for (i,j) in arcs_c2}
        M_c2=[i for i in range(1,int((len(nc2)+1)/2)+1)]
        arcs_c2_new, cost_c2_new = managing_pairs(arcs_c2,cost_c2)
        Ex1_c2_name,Ex1_c2_value, IS2 =IP_model_ind(nc2,arcs_c2_new,cost_c2_new,M_c2)
        print ('IS2',IS2)
        if str(IS2)=='None':
            IS2=0
        IS_c2.append(IS2)
        
        
        Node_c2=[]
        for i in range(len(Ex1_c2_name)):
            if Ex1_c2_name[i]!='__dummy':
                #print (Ex1_c2_name[i].split('_')[2].split('(')[1].split(',')[0],"=",Ex1_c2_value[i])
                Node_c2.append(int(Ex1_c2_name[i].split('_')[2].split('(')[1].split(',')[0]))
                BG_rec_c2.append(c2_rec_bg[int(Ex1_c2_name[i].split('_')[2].split('(')[1].split(',')[0])])
        Num_matched_c2.append(len(Node_c2))
        
        GR_BG_Rec_country1=GR_BG_Rec_country1+nc1_rec #Global Random Recipient country 1
        GR_BG_Rec_country2=GR_BG_Rec_country2+nc2_rec
        
        GR_BG_Don_country1=GR_BG_Don_country1+nc1_don
        GR_BG_Don_country2=GR_BG_Don_country2+nc2_don
        
        GR_BG_Rec_Total=GR_BG_Rec_country1+GR_BG_Rec_country2
        GR_BG_Don_Total=GR_BG_Don_country1+GR_BG_Don_country2
        
        nc1=[i for i in range(len(GR_BG_Rec_country1))] #Number of nodes in C1
        nc2=[i for i in range(len(GR_BG_Rec_country1),len(GR_BG_Rec_country1)+len(GR_BG_Rec_country2))] #Number of nodes in C2
        n=nc1+nc2 #Total number of nodes'
        
        arcs=compatible_edges(GR_BG_Rec_Total, GR_BG_Don_Total)
        c={(i,j):create_score(i,j) for (i,j) in arcs}
        
        M=[i for i in range(1,int(len(n)/2)+1)]
        
        #for individual rationality
        nc1_IR=[i for i in range(len(nc1_rec))] #Number of nodes in C1
        arcs_c1_IR=compatible_edges(nc1_rec,nc1_don)
        cost_c1_IR={(i,j):1 for (i,j) in arcs_c1_IR}
        M_c1_IR=[i for i in range(1,int(len(nc1_IR)/2))]
        arcs_c1_IR_new, cost_c1_IR_new = managing_pairs(arcs_c1_IR,cost_c1_IR)
        Ex1_c1_name_IR,Ex1_c1_value_IR, IS1_IR =IP_model_ind(nc1_IR,arcs_c1_IR_new,cost_c1_IR_new,M_c1_IR)
        print ('IS1_IR',IS1_IR)
        if str(IS1_IR)=='None':
            IS1_IR=0
        
        nc2_IR=[i for i in range(len(nc2_rec))] #Number of nodes in C2
        arcs_c2_IR=compatible_edges(nc2_rec,nc2_don)
        cost_c2_IR={(i,j):1 for (i,j) in arcs_c2_IR}
        M_c2_IR=[i for i in range(1,int((len(nc2_IR)+1)/2)+1)]
        arcs_c2_IR_new, cost_c2_IR_new = managing_pairs(arcs_c2_IR,cost_c2_IR)
        Ex1_c2_name_IR,Ex1_c2_value_IR, IS2_IR =IP_model_ind(nc2_IR,arcs_c2_IR_new,cost_c2_IR_new,M_c2_IR)
        print ('IS2_IR',IS2_IR)
        if str(IS2_IR)=='None':
            IS2_IR=0
        
        b1=int(IS1_IR)
        b2=int(IS2_IR)
        arcs_new, c_new = managing_pairs(arcs,c)
        
        Ex1_name, Ex1_value, GS_sol =IP_model_int(n,nc1,nc2,arcs_new,c_new,M,b1,b2)
        GS.append(GS_sol)
        
        Node=[]
        Matched_pair=[]
        for i in range(len(Ex1_name)):
            p=['__dummy','var1','var2']
            if Ex1_name[i] not in p:
                Node.append(int(Ex1_name[i].split('_')[2].split('(')[1].split(',')[0]))
                Matched_pair.append((int(Ex1_name[i].split('_')[2].split('(')[1].split(',')[0]),int(Ex1_name[i].split(',')[1].split('_')[1].split(')')[0])))
       
        C1_count=0
        C2_count=0
        print ('Matched pairs',Matched_pair)
        score_c1=[c[i,j] for (i,j) in Matched_pair if j in nc1]
        print ('score_c1',sum(score_c1))
        score_c2=[c[i,j] for (i,j) in Matched_pair if j not in nc1]
        print ('score_c2',sum(score_c2))
        for j in Node:
            if j in nc1:
                C1_count=C1_count+1
                GS_BG_rec_c1.append(GR_BG_Rec_Total[j])
            else:
                C2_count=C2_count+1
                GS_BG_rec_c2.append(GR_BG_Rec_Total[j])
        
        GS_c1.append(C1_count)
        GS_c2.append(C2_count)
        GS_c1_score.append(sum(score_c1))
        GS_c2_score.append(sum(score_c2))
        
        print ('c1_rec_bg,c1_don_bg',c1_rec_bg,c1_don_bg)
        print ('Node_c1',Node_c1)
        print ('c2_rec_bg,c2_don_bg',c2_rec_bg,c2_don_bg)
        print ('Node_c2',Node_c2)
        print ('Node matched',Node,len(Node))
        print ('n',n,len(n))
        
        #Program for deleting matched pairs and dropouts
        c1_rec_rem=[]
        c1_don_rem=[]
        index_c1=[]
        for i in range(len(c1_rec_bg)):
            if i not in Node_c1:
                random1=rng1.uniform(0,1)
                #print ('Random 1',random1)
                if random1 >=p1:
                    c1_rec_rem.append(c1_rec_bg[i])
                    c1_don_rem.append(c1_don_bg[i])
                    index_c1.append(i)
        print ('Remaining bg c1',c1_rec_rem)
        print ('Index C1',index_c1)
        c1_rec_bg=c1_rec_rem
        c1_don_bg=c1_don_rem
        
        c2_rec_rem=[]
        c2_don_rem=[]
        index_c2=[]
        for i in range(len(c2_rec_bg)):
            if i not in Node_c2:
                random2=rng1.uniform(0,1)
                #print ('Random 2',random2)
                if random2 >=p2:
                    c2_rec_rem.append(c2_rec_bg[i])
                    c2_don_rem.append(c2_don_bg[i])
                    index_c2.append(i)
        print ('Remaining bg c2',c2_rec_rem)
        print ('Index C1',index_c2)
        c2_rec_bg=c2_rec_rem
        c2_don_bg=c2_don_rem
        
        GR_c1_rec_rem=[]
        GR_c1_don_rem=[]
        GR_c2_rec_rem=[]
        GR_c2_don_rem=[]
        
        for i in range(len(GR_BG_Rec_Total)):
            #print ('i',i,'len(GR_BG_Rec_country1)',len(GR_BG_Rec_country1))
            if i<= len(GR_BG_Rec_country1)-1 and i not in Node:
                GS_random1=rng1.uniform(0,1)
                #print ('GS Random 1',GS_random1)
                if GS_random1 >=p1:
                    GR_c1_rec_rem.append(GR_BG_Rec_country1[i])
                    GR_c1_don_rem.append(GR_BG_Don_country1[i])
            elif i not in Node:
                GS_random2=rng1.uniform(0,1)
                #print ('GS Random 2',GS_random2)
                if GS_random2 >=p2:
                    GR_c2_rec_rem.append(GR_BG_Rec_country2[i-len(GR_BG_Rec_country1)])
                    GR_c2_don_rem.append(GR_BG_Don_country2[i-len(GR_BG_Rec_country1)])
        
        GR_BG_Rec_country1=GR_c1_rec_rem
        GR_BG_Don_country1=GR_c1_don_rem
        GR_BG_Rec_country2=GR_c2_rec_rem
        GR_BG_Don_country2=GR_c2_don_rem
        
        print ('Remaining Rec C1 in GS',GR_c1_rec_rem)
        print ('Remaining Rec C2 in GS',GR_c2_rec_rem)
        
    
    print ('IS_C1',IS_c1)
    print ('Individual Solution C1 over 20 rounds',sum(IS_c1),'mean',np.mean(IS_c1))
    print ('Individual Solution C2 over 20 rounds',sum(IS_c2),'mean',np.mean(IS_c2))
    print ('Sum of individual solution over 20 rounds',sum(IS_c1)+sum(IS_c2))
    print ('Global Solution over 20 rounds',sum(GS))
    print ('Share of C1 in GS',sum(GS_c1_score),'Share of C2 in GS',sum(GS_c2_score),'Total score',sum(GS_c1_score)+sum(GS_c2_score))
    print ('Number of patients matched for C1 in individual solution',sum(Num_matched_c1))
    print ('Number of patients matched for C2 in individual solution',sum(Num_matched_c2))
    print ('Number of patients matched for C1 in global solution',sum(GS_c1))
    print ('Number of patients matched for C2 in global solution',sum(GS_c2))
    
    print ('BG count of C1 in individual solution','O=',count_bg(BG_rec_c1,'O'),'A=',count_bg(BG_rec_c1,'A'),'B=',count_bg(BG_rec_c1,'B'),'AB=',count_bg(BG_rec_c1,'AB'))
    print ('BG count of C2 in individual solution','O=',count_bg(BG_rec_c2,'O'),'A=',count_bg(BG_rec_c2,'A'),'B=',count_bg(BG_rec_c2,'B'),'AB=',count_bg(BG_rec_c2,'AB'))
    print ('BG count of C1 in global solution','O=',count_bg(GS_BG_rec_c1,'O'),'A=',count_bg(GS_BG_rec_c1,'A'),'B=',count_bg(GS_BG_rec_c1,'B'),'AB=',count_bg(GS_BG_rec_c1,'AB'))
    print ('BG count of C2 in global solution','O=',count_bg(GS_BG_rec_c2,'O'),'A=',count_bg(GS_BG_rec_c2,'A'),'B=',count_bg(GS_BG_rec_c2,'B'),'AB=',count_bg(GS_BG_rec_c2,'AB'))     
    
    Avg_IS_C1.append(sum(Num_matched_c1))
    Avg_IS_C2.append(sum(Num_matched_c2))
    Avg_GS_C1.append(sum(GS_c1))
    Avg_GS_C2.append(sum(GS_c2))
    Avg_IS_C1_score.append(sum(IS_c1))
    Avg_IS_C2_score.append(sum(IS_c2))
    Avg_GS_C1_score.append(sum(GS_c1_score))
    Avg_GS_C2_score.append(sum(GS_c2_score))
    
    Avg_IS_C1_matched.append(Num_matched_c1)
    Avg_IS_C2_matched.append(Num_matched_c2)
    Avg_GS_C1_matched.append(GS_c1)
    Avg_GS_C2_matched.append(GS_c2)
    
    Avg_BG_IS_C1.append([count_bg(BG_rec_c1,'O'),count_bg(BG_rec_c1,'A'),count_bg(BG_rec_c1,'B'),count_bg(BG_rec_c1,'AB')])
    Avg_BG_IS_C2.append([count_bg(BG_rec_c2,'O'),count_bg(BG_rec_c2,'A'),count_bg(BG_rec_c2,'B'),count_bg(BG_rec_c2,'AB')])
    Avg_BG_GS_C1.append([count_bg(GS_BG_rec_c1,'O'),count_bg(GS_BG_rec_c1,'A'),count_bg(GS_BG_rec_c1,'B'),count_bg(GS_BG_rec_c1,'AB')])
    Avg_BG_GS_C2.append([count_bg(GS_BG_rec_c2,'O'),count_bg(GS_BG_rec_c2,'A'),count_bg(GS_BG_rec_c2,'B'),count_bg(GS_BG_rec_c2,'AB')])
# =============================================================================
#     print ('BG of C1 in individual solution',BG_rec_c1,len(BG_rec_c1))
#     print ('BG of C2 in individual solution',BG_rec_c2,len(BG_rec_c2))
#     print ('BG of C1 in global solution',GS_BG_rec_c1,len(GS_BG_rec_c1))
#     print ('BG of C2 in global solution',GS_BG_rec_c2,len(GS_BG_rec_c2))
#     
#     print ('BG count of C1 in individual solution','O=',count_bg(BG_rec_c1,'O'),'A=',count_bg(BG_rec_c1,'A'),'B=',count_bg(BG_rec_c1,'B'),'AB=',count_bg(BG_rec_c1,'AB'))
#     print ('BG count of C2 in individual solution','O=',count_bg(BG_rec_c2,'O'),'A=',count_bg(BG_rec_c2,'A'),'B=',count_bg(BG_rec_c2,'B'),'AB=',count_bg(BG_rec_c2,'AB'))
#     print ('BG count of C1 in global solution','O=',count_bg(GS_BG_rec_c1,'O'),'A=',count_bg(GS_BG_rec_c1,'A'),'B=',count_bg(GS_BG_rec_c1,'B'),'AB=',count_bg(GS_BG_rec_c1,'AB'))
#     print ('BG count of C2 in global solution','O=',count_bg(GS_BG_rec_c2,'O'),'A=',count_bg(GS_BG_rec_c2,'A'),'B=',count_bg(GS_BG_rec_c2,'B'),'AB=',count_bg(GS_BG_rec_c2,'AB'))
#     
#         
#     print ('Total time in sec', time.clock() - start_time)
#     
#     with open("Output.txt", "w") as text_file:
#         print ('Solution for C1-U(5-10), C2-U(5,10), prob_c1=0.2,prob_c2=0.2,rounds=20')
#         print('Individual Solution C1 over 20 rounds', sum(IS_c1),'mean', np.mean(IS_c1), file=text_file)
#         print ('Individual Solution C2 over 20 rounds',sum(IS_c2),'mean',np.mean(IS_c2),file=text_file)
#         print ('Sum of individual solution over 20 rounds',sum(IS_c1)+sum(IS_c2),file=text_file)
#         print ('Global Solution over 20 rounds',sum(GS),file=text_file)
#         print ('Number of patients matched for C1 in individual solution',sum(Num_matched_c1), file=text_file)
#         print ('Number of patients matched for C2 in individual solution',sum(Num_matched_c2), file=text_file)
#         print ('Number of patients matched for C1 in global solution',sum(GS_c1), file=text_file)
#         print ('Number of patients matched for C2 in global solution',sum(GS_c2), file=text_file)
#         print ('BG count of C1 in individual solution','O=',count_bg(BG_rec_c1,'O'),'A=',count_bg(BG_rec_c1,'A'),'B=',count_bg(BG_rec_c1,'B'),'AB=',count_bg(BG_rec_c1,'AB'), file=text_file)
#         print ('BG count of C2 in individual solution','O=',count_bg(BG_rec_c2,'O'),'A=',count_bg(BG_rec_c2,'A'),'B=',count_bg(BG_rec_c2,'B'),'AB=',count_bg(BG_rec_c2,'AB'), file=text_file)
#         print ('BG count of C1 in global solution','O=',count_bg(GS_BG_rec_c1,'O'),'A=',count_bg(GS_BG_rec_c1,'A'),'B=',count_bg(GS_BG_rec_c1,'B'),'AB=',count_bg(GS_BG_rec_c1,'AB'), file=text_file)
#         print ('BG count of C2 in global solution','O=',count_bg(GS_BG_rec_c2,'O'),'A=',count_bg(GS_BG_rec_c2,'A'),'B=',count_bg(GS_BG_rec_c2,'B'),'AB=',count_bg(GS_BG_rec_c2,'AB'), file=text_file)
#     
#         xpos=[i for i in range(z)]
#         plt.plot(xpos,Num_matched_c1,'r--',label='Individual solution C1')
#         plt.plot(xpos,Num_matched_c2,'b--',label='Individual solution C2')
#         plt.plot(xpos,GS_c1,'r',label='Global solution C1')
#         plt.plot(xpos,GS_c2,'b',label='Global solution C2')
#         plt.xlabel('Number of Rounds')
#         plt.ylabel('Number of recipient matched')
#         plt.title('Comparison of Individual solution vs Global solution for 2 countries')
#         plt.legend()
#         plt.savefig('Multi-registry.png')
# 
# =============================================================================
print ('Avg #matches for C1 in Ind solution',np.mean(Avg_IS_C1),'Avg Score',np.mean(Avg_IS_C1_score))
print ('Avg #matches for C2 in Ind solution',np.mean(Avg_IS_C2),'Avg Score',np.mean(Avg_IS_C2_score))
print ('Avg #matches for C1 in Global solution',np.mean(Avg_GS_C1),'Avg Score',np.mean(Avg_GS_C1_score))
print ('Avg #matches for C2 in Global solution',np.mean(Avg_GS_C2),'Avg Score',np.mean(Avg_GS_C2_score))
       
#print ('Avg_IS_C1_matched',Avg_IS_C1_matched)
#print ('Avg_IS_C2_matched',Avg_IS_C2_matched)
#print ('Avg_GS_C1_matched',Avg_GS_C1_matched)
#print ('Avg_GS_C2_matched',Avg_GS_C2_matched)

#print ('Avg_BG_IS_C1',Avg_BG_IS_C1)
#print ('Avg_BG_IS_C2',Avg_BG_IS_C2)
#print ('Avg_BG_GS_C1',Avg_BG_GS_C1)
#print ('Avg_BG_GS_C2',Avg_BG_GS_C2)

Avg_IS_c1=[] #Average number of transplants for each round 
Avg_IS_c2=[]
Avg_GS_c1=[]
Avg_GS_c2=[]
      
for j in range(runs):
    IS_c1=[]
    IS_c2=[]
    GS_c1=[]
    GS_c2=[]
    for i in range(rep):
        IS_c1.append(Avg_IS_C1_matched[i][j])
        IS_c2.append(Avg_IS_C2_matched[i][j])
        GS_c1.append(Avg_GS_C1_matched[i][j])
        GS_c2.append(Avg_GS_C2_matched[i][j])
    Avg_IS_c1.append(np.mean(IS_c1))
    Avg_IS_c2.append(np.mean(IS_c2))
    Avg_GS_c1.append(np.mean(GS_c1))
    Avg_GS_c2.append(np.mean(GS_c2))
print ('Avg_IS_c1',Avg_IS_c1)
print ('Avg_IS_c2',Avg_IS_c2)
print ('Avg_GS_c1',Avg_GS_c1)
print ('Avg_GS_c2',Avg_GS_c2)
       
BG_IS_C1=[]
BG_IS_C2=[]
BG_GS_C1=[]
BG_GS_C2=[]
for j in range(4):
    count_IS_C1=[]
    count_IS_C2=[]
    count_GS_C1=[]
    count_GS_C2=[]
    for i in range(rep):
        count_IS_C1.append(Avg_BG_IS_C1[i][j])
        count_IS_C2.append(Avg_BG_IS_C2[i][j])
        count_GS_C1.append(Avg_BG_GS_C1[i][j])
        count_GS_C2.append(Avg_BG_GS_C2[i][j])
    BG_IS_C1.append(np.mean(count_IS_C1))
    BG_IS_C2.append(np.mean(count_IS_C2))
    BG_GS_C1.append(np.mean(count_GS_C1))
    BG_GS_C2.append(np.mean(count_GS_C2))

print ('BG_IS_C1',BG_IS_C1,'O, A, B, AB')
print ('BG_IS_C2',BG_IS_C2,'O, A, B, AB')
print ('BG_GS_C1',BG_GS_C1,'O, A, B, AB')
print ('BG_GS_C2',BG_GS_C2,'O, A, B, AB')

# =============================================================================
# xpos=[i for i in range(runs)]
# plt.plot(xpos,Avg_IS_c1,'r--',label='Individual solution C1')
# plt.plot(xpos,Avg_IS_c2,'b--',label='Individual solution C2')
# plt.plot(xpos,Avg_GS_c1,'r',label='Global solution C1')
# plt.plot(xpos,Avg_GS_c2,'b',label='Global solution C2')
# plt.xlabel('Number of Rounds')
# plt.ylabel('Number of recipient matched')
# plt.title('Comparison of Individual solution vs Global solution for 2 countries')
# plt.legend()
# plt.savefig('Multi-registry.png')
# 
# =============================================================================
with open("Output.txt", "w") as text_file:
    print ('Solution for C1-U(5-10), C2-U(5,10), prob_c1=0.2,prob_c2=0.2,','rounds=',runs,'replications',rep)
    print ('Avg #matches for C1 in Ind solution',np.mean(Avg_IS_C1),'Avg Score',np.mean(Avg_IS_C1_score))
    print ('Avg #matches for C2 in Ind solution',np.mean(Avg_IS_C2),'Avg Score',np.mean(Avg_IS_C2_score))
    print ('Avg #matches for C1 in Global solution',np.mean(Avg_GS_C1),'Avg Score',np.mean(Avg_GS_C1_score))
    print ('Avg #matches for C2 in Global solution',np.mean(Avg_GS_C2),'Avg Score',np.mean(Avg_GS_C2_score))
    print ('BG_IS_C1',BG_IS_C1,'O, A, B, AB')
    print ('BG_IS_C2',BG_IS_C2,'O, A, B, AB')
    print ('BG_GS_C1',BG_GS_C1,'O, A, B, AB')
    print ('BG_GS_C2',BG_GS_C2,'O, A, B, AB') 
    
    xpos=[i for i in range(runs)]
    plt.plot(xpos,Avg_IS_c1,'r--',label='Individual solution C1')
    plt.plot(xpos,Avg_IS_c2,'b--',label='Individual solution C2')
    plt.plot(xpos,Avg_GS_c1,'r',label='Global solution C1')
    plt.plot(xpos,Avg_GS_c2,'b',label='Global solution C2')
    plt.xlabel('Number of Rounds \n c1-u(5,10)')
    plt.ylabel('Number of recipient matched')
    #plt.title('Comparison of Individual solution vs Global solution for 2 countries')
    #plt.figure(facecolor="white")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=1.)
    plt.savefig('Multi-registry.png')

