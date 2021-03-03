# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 21:11:22 2020

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
import pulp as lp

rng1 = np.random.RandomState(11)



def IP_model_DDIC(n,P,PWL,DD,WL,arcs,c):
    
    M=[i for i in range(1,int(len(n)/2)+1)]
    #print ('arcs',arcs,'c',c)
    prob1= LpProblem("Int_model_for_KEP", LpMaximize)
    x = LpVariable.dicts('Route var',(arcs,M), cat='Binary')
    
    prob1 += lpSum([x[i][j]*c[i] for i in arcs for j in M])
    
    for k in P:
        for j in M:
            prob1 += lpSum([x[i][j] for i in arcs if i[0]==k]) == lpSum([x[i][j] for i in arcs if i[1]==k]) 
            
    for k in PWL:
        for j in M:
            prob1 += lpSum([x[i][j] for i in arcs if i[0]==k]) <= lpSum([x[i][j] for i in arcs if i[1]==k])
    
    for k in n:
        prob1 += lpSum([x[i][j] for i in arcs for j in M if i[0]==k]) <=1
    
    for k in n:
        prob1 += lpSum([x[i][j] for i in arcs for j in M if i[1]==k]) <=1
        
    for l in M:
        prob1 += lpSum([x[i][l] for i in arcs]) <=3
            
    prob1.writeLP('Matching.LP')
    
    prob1.solve()
    print ('Status',LpStatus[prob1.status])
    print ('Objective DDIC',lp.value(prob1.objective))
    #print ('x',prob1.variables())
    
    z_name=[]
    z_value=[]
    count=0
    for i in prob1.variables():
        if i.value()!=0:
            z_name.append(i.name)
            z_value.append(i.value())
            count=count+1
            print (i,i.value())
    
    return z_name,z_value,lp.value(prob1.objective)

def IP_model_ind(n,arcs,c):
        
    #print ('arcs',arcs,'c',c)
    M=[i for i in range(1,int(len(n)/2)+1)]
    prob1= LpProblem("Int_model_for_KEP", LpMaximize)
    x = LpVariable.dicts('Route var',(arcs,M),0,1, cat='Integer')
      
    prob1 += lpSum([x[i][j]*c[i] for i in arcs for j in M])

    for k in n:
        for j in M:
            prob1 += lpSum([x[i][j] for i in arcs if i[0]==k]) - lpSum([x[i][j] for i in arcs if i[1]==k]) ==0
    
    for k in n:
        prob1 += lpSum([x[i][j] for i in arcs for j in M if i[0]==k]) <=1
        
    for l in M:
        prob1 += lpSum([x[i][l] for i in arcs]) <=3
            
    prob1.writeLP('Matching.LP')
    
    prob1.solve()
    print ('Status GS',LpStatus[prob1.status])
    print ('Objective Ind',lp.value(prob1.objective))
    #print ('x',prob1.variables())
    
    z_name=[]
    z_value=[]
    count=0
    for i in prob1.variables():
        if i.value()!=0:
            z_name.append(i.name)
            z_value.append(i.value())
            count=count+1
            print (i,i.value())
    
    return z_name,z_value,lp.value(prob1.objective)


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

#Generating BG for DD
def generating_bg_DD(n):
	x=[0.37,0.69,0.92,1] #BG probability - O,B,A,AB
	generated_DD=[]
	for i in range(n):
		random_num=rng1.uniform(0,1)
		if random_num <=x[0]:
			generated_DD.append('O')
		elif random_num >= x[0] and random_num <= x[1]:
			generated_DD.append('B')
		elif random_num >= x[1] and random_num <= x[2]:
			generated_DD.append('A')
		elif random_num >= x[2] and random_num <= 1:
			generated_DD.append('AB')
			
	return generated_DD 

def compatible_edges(recipient_bg, donor_bg):
    Edges=[]
    for i in range(0,len(donor_bg)):
        for j in range(0,len(recipient_bg)):
            if ((donor_bg[i]==recipient_bg[j] and i!=j or donor_bg[i]=='O' or recipient_bg[j]=='AB' or recipient_bg[j]=='WL') and (donor_bg[i]!="Matched" and recipient_bg[i]!="Matched" and donor_bg[j]!="Matched" and recipient_bg[j]!="Matched")):
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
    age_score=0
    if donor_age - rec_age <= -20 or donor_age - rec_age >= 20 :
        age_score=15
    if (donor_age - rec_age >= -20 and donor_age - rec_age <= -10) or (donor_age - rec_age >= 10 and donor_age - rec_age <= 20) :
        age_score=30
    if (donor_age - rec_age >= -10 and donor_age - rec_age <= 0) or (donor_age - rec_age >= 0 and donor_age - rec_age <= 10) :
        age_score=50
        
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

def managing_pairs(arcs,cost):
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

#Example of DDIC

n=[1,2,3,4,5,6,7]
DD=[1,2]
P=[3,4]
PWL=[5]
WL=[6,7]

arcs=[(1,3),(2,5),(3,4),(3,6),(3,7),(4,6),(4,7),(5,6),(5,7)]
cost={(i,j):1 for (i,j) in arcs}

name, value, objective = IP_model_DDIC(n,P,PWL,DD,WL,arcs,cost)

print ('Objective value',objective)
print ('Variable name',name,'value =',value)

#Program for multiple replications

Avg_DDIC_sol=[]
Avg_Ind_sol=[]
Avg_DDIC_objective=[]
Avg_Ind_objective=[]
Avg_DDIC_dropout=[]
Avg_Ind_dropout=[]


Avg_DDIC_matched=[]
Avg_Ind_matched=[]
Avg_DDIC_DO=[]
Avg_Ind_DO=[]
Avg_DDIC_wait_time=[]
Avg_Ind_wait_time=[]

Avg_DDIC_BG_waiting_time=[]
Avg_Ind_BG_waiting_time=[]

Avg_BG_DDIC=[]
Avg_BG_Ind=[]
Avg_Just_DD_sol=[]
Avg_Just_Ind_sol=[]
Effective_DO_prob_DDIC=[] #To calculate effective dropout probability in DDIC
Effective_DO_prob_Ind=[] 

Replications=20
Runs=12

for b in range(Replications):
    
    Num_pair_generated=[] #Number of pairs generated in each round
    Num_DD_generated=[] #Number of DD generated in each round
    
    Donor_BG_total=[] #List of donor's BG for pairs
    Recipient_BG_total=[] #List of recipient's BG for pairs
    
    DDIC_Donor_BG_pair_total=[] #List of donor's BG for pairs in DDIC allocation
    DDIC_Recipient_BG_pair_total=[] #List of recipient's BG for pairs in DDIC allocation
    
    DD_BG_total=[] #List of DD's BG
    
    Num_pair_matched=[] #Number of pairs matched in each round
    Num_WL_matched=[] #Number of wait list patients matched in each round
    
    Additional_gain=[]
    Swap_solution=[]
    DDIC_solution=[]
    
    
    Pair_waiting_time=0 # Waiting time in swap allocation
    Dropout_pair=[] #Number of dropout pairs in swap allocation
    DD_Pair_waiting_time=0
    DD_Dropout_pair=[]
    DDIC_wait_time=[]
    Ind_wait_time=[]
    
    WT_BG_O_Swap=0 #Waiting time for O recipient in swap registry
    WT_BG_A_Swap=0
    WT_BG_B_Swap=0
    WT_BG_AB_Swap=0
    
    WT_BG_O_DDIC=0 #Waiting time for O recipient in DDIC
    WT_BG_A_DDIC=0
    WT_BG_B_DDIC=0
    WT_BG_AB_DDIC=0
    
    DO_BG_O_Swap=0 # Total number of dropouts for O recipients in swap registry
    DO_BG_A_Swap=0
    DO_BG_B_Swap=0
    DO_BG_AB_Swap=0
    
    DO_BG_O_DDIC=0 #Total number of dropouts for O recipient in DDIC
    DO_BG_A_DDIC=0
    DO_BG_B_DDIC=0
    DO_BG_AB_DDIC=0
    
    count_O=0 #Count total number of O recipient generated in swap registry
    count_A=0
    count_B=0
    count_AB=0
    
    p=0.4 #probability of dropout for a pair
    z=Runs
    
    donor_bg_P_total=[]
    donor_bg_PWL_total=[]
    recipient_bg_P_total=[]
    recipient_bg_PWL_total=[]
    Num_matched=[]
    Num_matched_ind=[]
    
    DDIC_BG_waiting_time=[]
    Ind_BG_waiting_time=[]
    
    Recipient_BG_total_ind=[]
    Donor_BG_total_ind=[]
    
    BG_DDIC=[]
    BG_Ind=[]
    Just_DD_sol=[]
    Just_Ind_sol=[]
    
    Eff_DO_Prob_DDIC=[]
    Eff_DO_Prob_Ind=[]
    
    Objective_DDIC=[]
    Objective_Ind=[]
    
    for a in range(z):
        print ('Round',a)
        random_DD=rng1.randint(1,5) #Number of DD generated in current round
        random_P=rng1.randint(1,5) #Number of Pair generated in current round
        random_PWL=rng1.randint(1,5)
        dd_blood_group=generating_bg_DD(random_DD) #BG generation of DD
        print ('DD BG',dd_blood_group)
        recipient_bg_P, donor_bg_P = generating_bg_ASTRA(random_P) #BG generation of Pairs
        print ('Donor BG Pair',donor_bg_P,'Rec BG Pair',recipient_bg_P)
        recipient_bg_PWL, donor_bg_PWL = generating_bg_ASTRA(random_PWL) #BG generation of Pairs
        print ('Donor BG PWL',donor_bg_PWL,'Rec BG PWL',recipient_bg_PWL)
        
        WL_patients=[]
        for i in range(random_DD):
            WL_patients.append('WL')
        #program for DDIC    
        recipient_bg_P_total=recipient_bg_P_total+recipient_bg_P
        recipient_bg_PWL_total=recipient_bg_PWL_total+recipient_bg_PWL
        donor_bg_P_total=donor_bg_P_total+donor_bg_P
        donor_bg_PWL_total=donor_bg_PWL_total+donor_bg_PWL
        
        Recipient_BG_total=recipient_bg_P_total+recipient_bg_PWL_total+WL_patients
        Donor_BG_total=donor_bg_P_total+donor_bg_PWL_total+dd_blood_group
        
        arcs=compatible_edges(Recipient_BG_total,Donor_BG_total)
        cost={(i,j):create_score(i,j) for (i,j) in arcs}
        arcs_new, cost_new = managing_pairs(arcs,cost)
        
        n=[i for i in range(len(Recipient_BG_total))]
        P=[i for i in range(len(recipient_bg_P_total))]
        PWL=[i for i in range(len(recipient_bg_P_total),len(recipient_bg_P_total)+len(recipient_bg_PWL_total))]
        DD=[i for i in range(len(donor_bg_P_total)+len(donor_bg_PWL_total),len(donor_bg_P_total)+len(donor_bg_PWL_total)+len(dd_blood_group))]
        WL=[i for i in range(len(recipient_bg_P_total)+len(recipient_bg_PWL_total),len(recipient_bg_P_total)+len(recipient_bg_PWL_total)+len(WL_patients))]
        
        print ('n',n,'P',P,'PWL',PWL,'DD',DD,'WL',WL)
        name, value, objective = IP_model_DDIC(n,P,PWL,DD,WL,arcs_new,cost_new)
        
        print ('objective DDIC',objective)
        #print ('Name',name)
        if str(objective)=='None':
            objective=0
        Node=[]
        for i in range(len(name)):
            if name[i]!='__dummy':
                print (name[i].split('_')[2].split('(')[1].split(',')[0],"=",value[i])
                Node.append(int(name[i].split('_')[2].split('(')[1].split(',')[0]))
    
        Num_matched.append(len(Node))
        Objective_DDIC.append(objective)
        
        k=0 #Total number of P+PWL pairs matched in each round
        for j in Node:
            if j in P+PWL:
                BG_DDIC.append(Recipient_BG_total[j])
                k=k+1
            else:
                Just_DD_sol.append(Recipient_BG_total[j])
                
        #print ('BG_DDIC',BG_DDIC)
        #print ('DD_sol',Just_DD_sol)
        #print ('BG DDIC, P,PWL, DD, n',k,len(P),len(PWL),len(DD),len(n))
        Eff_DO_Prob_DDIC.append(1-(k/(len(n)-len(DD)))) #proportion of dropouts in each round
        
        #Deleting dropout pairs
        rec_rem=[]
        don_rem=[]
        index=[]
        count_dropout_P=0
        Wait_time_P=0
        for i in range(len(recipient_bg_P_total)):
            if i not in Node:
                random1=rng1.uniform(0,1)
                #print ('Random 1',random1)
                if random1 >=p:
                    rec_rem.append(recipient_bg_P_total[i])
                    don_rem.append(donor_bg_P_total[i])
                    index.append(i)
                    Wait_time_P=Wait_time_P+3
                else:
                    count_dropout_P=count_dropout_P+1
        print ('Remaining bg Pair',rec_rem)
        print ('Index P',index)
        recipient_bg_P_total=rec_rem
        donor_bg_P_total=don_rem
        DDIC_BG_waiting_time=DDIC_BG_waiting_time+rec_rem
        
        rec_rem=[]
        don_rem=[]
        index=[]
        count_dropout_PWL=0
        Wait_time_PWL=0
        for i in range(len(recipient_bg_P_total),len(recipient_bg_P_total)+len(recipient_bg_PWL_total)):
            if i not in Node:
                random1=rng1.uniform(0,1)
                #print ('Random 1',random1)
                if random1 >=p:
                    rec_rem.append(recipient_bg_PWL_total[i-len(recipient_bg_P_total)])
                    don_rem.append(donor_bg_PWL_total[i-len(donor_bg_P_total)])
                    index.append(i)
                    Wait_time_PWL=Wait_time_PWL+3
                else:
                    count_dropout_PWL=count_dropout_PWL+1
        print ('Remaining bg Pair WL',rec_rem)
        print ('Index PWL',index)
        recipient_bg_PWL_total=rec_rem
        donor_bg_PWL_total=don_rem
        DDIC_BG_waiting_time=DDIC_BG_waiting_time+rec_rem
        
        DD_Dropout_pair.append(count_dropout_P+count_dropout_PWL)
        DDIC_wait_time.append(Wait_time_PWL+Wait_time_P)
        
        #program for Independent functioning of registry
        Recipient_BG_total_ind=Recipient_BG_total_ind+recipient_bg_P+recipient_bg_PWL
        Donor_BG_total_ind=Donor_BG_total_ind+donor_bg_P+donor_bg_PWL
        
        n=[i for i in range(len(Recipient_BG_total_ind))]
        arcs=compatible_edges(Recipient_BG_total_ind,Donor_BG_total_ind)
        cost={(i,j):create_score(i,j) for (i,j) in arcs}
        arcs_new, cost_new = managing_pairs(arcs,cost)
        name, value, objective = IP_model_ind(n,arcs_new,cost_new)
        
        print ('objective ind',objective)
        #print ('Name ind',name)
        if str(objective)=='None':
            objective=0
        Node=[]
        for i in range(len(name)):
            if name[i]!='__dummy':
                #print('hello')
                print (name[i].split('_')[2].split('(')[1].split(',')[0],"=",value[i])
                Node.append(int(name[i].split('_')[2].split('(')[1].split(',')[0]))
                BG_Ind.append(Recipient_BG_total_ind[int(name[i].split('_')[2].split('(')[1].split(',')[0])])
        
        Num_matched_ind.append(len(Node)+random_DD)
        Objective_Ind.append(objective+75*random_DD)
        Just_Ind_sol.append(len(Node))
        
        print ('len(node), len(n)',len(Node),len(n))
        Eff_DO_Prob_Ind.append(1-(len(Node)/len(n)))
        
        #Deleting dropout pairs
        rec_rem=[]
        don_rem=[]
        index=[]
        count_dropout=0
        Wait_time_ind=0
        for i in range(len(Recipient_BG_total_ind)):
            if i not in Node:
                random1=rng1.uniform(0,1)
                #print ('Random 1',random1)
                if random1 >=p:
                    rec_rem.append(Recipient_BG_total_ind[i])
                    don_rem.append(Donor_BG_total_ind[i])
                    index.append(i)
                    Wait_time_ind=Wait_time_ind+3
                else:
                    count_dropout=count_dropout+1
        print ('Remaining Pair in Ind functioning registry',rec_rem,don_rem)
        print ('Index ind',index)
        Recipient_BG_total_ind=rec_rem
        Donor_BG_total_ind=don_rem
        Ind_BG_waiting_time=Ind_BG_waiting_time+rec_rem
        Dropout_pair.append(count_dropout)
        Ind_wait_time.append(Wait_time_ind)
        
    print ('Total number of recipient matched over 12 rounds in DDIC',sum(Num_matched),'Average',np.mean(Num_matched)) 
    print ('Total number of recipient matched over 12 rounds in Ind functioning',sum(Num_matched_ind),'Average',np.mean(Num_matched_ind)) 
    
    print ('Total score in DDIC over 12 rounds',sum(Objective_DDIC),'mean',np.mean(Objective_DDIC))
    print ('Total score in Ind functioning registry over 12 rounds',sum(Objective_Ind),'mean',np.mean(Objective_Ind))
    
    print ('Total number of dropouts in DDIC over 12 rounds',sum(DD_Dropout_pair),'mean',np.mean(DD_Dropout_pair))
    print ('Total number of dropouts in Ind functioning registry over 12 rounds',sum(Dropout_pair),'mean',np.mean(Dropout_pair))
    
    Avg_DDIC_sol.append(sum(Num_matched))
    Avg_Ind_sol.append(sum(Num_matched_ind))
    Avg_DDIC_objective.append(sum(Objective_DDIC))
    Avg_Ind_objective.append(sum(Objective_Ind))
    Avg_DDIC_dropout.append(sum(DD_Dropout_pair))
    Avg_Ind_dropout.append(sum(Dropout_pair))
    
    Avg_DDIC_matched.append(Num_matched)
    Avg_Ind_matched.append(Num_matched_ind)
    Avg_DDIC_DO.append(DD_Dropout_pair)
    Avg_Ind_DO.append(Dropout_pair)
    print ('DDIC_wait_time',DDIC_wait_time)
    Avg_DDIC_wait_time.append(DDIC_wait_time)
    Avg_Ind_wait_time.append(Ind_wait_time)
    Avg_DDIC_BG_waiting_time.append([count_bg(DDIC_BG_waiting_time,'O'),count_bg(DDIC_BG_waiting_time,'A'),count_bg(DDIC_BG_waiting_time,'B'),count_bg(DDIC_BG_waiting_time,'AB')])
    Avg_Ind_BG_waiting_time.append([count_bg(Ind_BG_waiting_time,'O'),count_bg(Ind_BG_waiting_time,'A'),count_bg(Ind_BG_waiting_time,'B'),count_bg(Ind_BG_waiting_time,'AB')])
    print ('DDIC waiting time for different BG, O',count_bg(DDIC_BG_waiting_time,'O'),'A',count_bg(DDIC_BG_waiting_time,'A'),'B',count_bg(DDIC_BG_waiting_time,'B'),'AB',count_bg(DDIC_BG_waiting_time,'AB'))
    print ('Ind waiting time for different BG, O',count_bg(Ind_BG_waiting_time,'O'),'A',count_bg(Ind_BG_waiting_time,'A'),'B',count_bg(Ind_BG_waiting_time,'B'),'AB',count_bg(Ind_BG_waiting_time,'AB'))
    
    print ('BG count DDIC','O=',count_bg(BG_DDIC,'O'),'A=',count_bg(BG_DDIC,'A'),'B=',count_bg(BG_DDIC,'B'),'AB=',count_bg(BG_DDIC,'AB'))
    print ('DD BG count=', count_bg(Just_DD_sol,'WL'))
    print ('BG count Ind','O=',count_bg(BG_Ind,'O'),'A=',count_bg(BG_Ind,'A'),'B=',count_bg(BG_Ind,'B'),'AB=',count_bg(BG_Ind,'AB'))
    Avg_BG_DDIC.append([count_bg(BG_DDIC,'O'),count_bg(BG_DDIC,'A'),count_bg(BG_DDIC,'B'),count_bg(BG_DDIC,'AB')])
    Avg_BG_Ind.append([count_bg(BG_Ind,'O'),count_bg(BG_Ind,'A'),count_bg(BG_Ind,'B'),count_bg(BG_Ind,'AB')])
    Avg_Just_DD_sol.append(count_bg(Just_DD_sol,'WL'))
    Avg_Just_Ind_sol.append(sum(Just_Ind_sol))
    
    print ('Effective DO probability in DDIC over 12 rounds',Eff_DO_Prob_DDIC)
    Effective_DO_prob_DDIC.append(np.mean(Eff_DO_Prob_DDIC))
    Effective_DO_prob_Ind.append(np.mean(Eff_DO_Prob_Ind))
    # =============================================================================
    # xpos=[i for i in range(z)]
    # plt.plot(xpos,Num_matched,'r',label='DDIC Solution')
    # plt.plot(xpos,Num_matched_ind,'b',label='Independent functioning registries')
    # plt.plot(xpos,DD_Dropout_pair,'r--',label='DDIC dropouts')
    # plt.plot(xpos,Dropout_pair,'b--',label='Independent registry dropouts')
    # plt.xlabel('Number of Rounds')
    # plt.ylabel('Number of recipient matched')
    # plt.title('Comparison of DDIC solution vs Independent solution')
    # plt.legend()
    # plt.show()
    # =============================================================================

print ('Average DDIC sol',np.mean(Avg_DDIC_sol))
print ('Average Ind sol',np.mean(Avg_Ind_sol))
print ('Average DDIC objective',np.mean(Avg_DDIC_objective))
print ('Average Ind objective',np.mean(Avg_Ind_objective))
print ('Average DDIC dropouts',np.mean(Avg_DDIC_dropout))
print ('Average Ind dropouts',np.mean(Avg_Ind_dropout))
    
    
DDIC_sol_per_round=[] #Average number of transplants for each round 
Ind_sol_per_round=[]
DDIC_DO_per_round=[]
Ind_DO_per_round=[]
DDIC_wait_time_per_round=[]
Ind_wait_time_per_round=[]
      
for j in range(Runs):
    IS_c1=[]
    IS_c2=[]
    GS_c1=[]
    GS_c2=[]
    DDIC_wait=[]
    Ind_wait=[]
    for i in range(Replications):
        IS_c1.append(Avg_DDIC_matched[i][j])
        IS_c2.append(Avg_Ind_matched[i][j])
        GS_c1.append(Avg_DDIC_DO[i][j])
        GS_c2.append(Avg_Ind_DO[i][j])
        DDIC_wait.append(Avg_DDIC_wait_time[i][j])
        Ind_wait.append(Avg_Ind_wait_time[i][j])
        
    DDIC_sol_per_round.append(np.mean(IS_c1))
    Ind_sol_per_round.append(np.mean(IS_c2))
    DDIC_DO_per_round.append(np.mean(GS_c1))
    Ind_DO_per_round.append(np.mean(GS_c2))
    DDIC_wait_time_per_round.append(np.mean(DDIC_wait))
    Ind_wait_time_per_round.append(np.mean(Ind_wait))
    
#print ('Avg DDIC sol per round',DDIC_sol_per_round)
#print ('Avg Ind sol per round',Ind_sol_per_round)
print ('Avg DDIC DO per round',np.mean(DDIC_DO_per_round))
print ('Avg Ind DO per round',np.mean(Ind_DO_per_round))  
#print ('Avg DDIC waiting time per round',DDIC_wait_time_per_round)
#print ('Avg Ind waiting time per round',Ind_wait_time_per_round)
print ('Total waiting time DDIC in months', sum(DDIC_wait_time_per_round))
print ('Total waiting time Ind sol in months',sum(Ind_wait_time_per_round))

DDIC_BG_wait=[]
Ind_BG_wait=[]
DDIC_BG_matched=[]
Ind_BG_matched=[]
for j in range(4):
    count_DDIC=[]
    count_Ind=[]
    matched_DDIC=[]
    matched_Ind=[]
    for i in range(Replications):
        count_DDIC.append(Avg_DDIC_BG_waiting_time[i][j])
        count_Ind.append(Avg_Ind_BG_waiting_time[i][j])
        matched_DDIC.append(Avg_BG_DDIC[i][j])
        matched_Ind.append(Avg_BG_Ind[i][j])
        
    DDIC_BG_wait.append(np.mean(count_DDIC))
    Ind_BG_wait.append(np.mean(count_Ind))
    DDIC_BG_matched.append(np.mean(matched_DDIC))
    Ind_BG_matched.append(np.mean(matched_Ind))
    

print ('BG_waiting_DDIC',[i for i in DDIC_BG_wait],'O, A, B, AB')
print ('BG_waiting_Ind',[i for i in Ind_BG_wait],'O, A, B, AB')
print ('BG_matched_DDIC',[i for i in DDIC_BG_matched],'O, A, B, AB', 'total',sum(DDIC_BG_matched))
print ('BG_matched_Ind',[i for i in Ind_BG_matched],'O, A, B, AB', 'total',sum(Ind_BG_matched))
print ('Avg DD sol only',np.mean(Avg_Just_DD_sol))
print ('Avg Ind sol only',np.mean(Avg_Just_Ind_sol))

print ('Results to be written in manuscript')

print ('Average DDIC sol=',np.mean(Avg_DDIC_sol))
print ('Average Ind sol=',np.mean(Avg_Ind_sol))
print ('%gain in DDIC over Ind=',round((np.mean(Avg_DDIC_sol)-np.mean(Avg_Ind_sol))*100/np.mean(Avg_Ind_sol),1),'%')

Edge_score_DDIC=np.mean(Avg_DDIC_objective)/np.mean(Avg_DDIC_sol)
Edge_score_Ind=np.mean(Avg_Ind_objective)/np.mean(Avg_Ind_sol)
print ('Avg edge score DDIC=',round(Edge_score_DDIC,1))
print ('Avg edge score Ind=',round(Edge_score_Ind,1))
print ('%gain in avg edge score in DDIC over Ind=',round((Edge_score_DDIC-Edge_score_Ind)*100/Edge_score_Ind,1),'%')
a=np.mean(Avg_DDIC_sol)-np.mean(Avg_Just_DD_sol)
print ('Avg waiting time per recipient in DDIC',sum(DDIC_wait_time_per_round)/a)
print ('Avg waiting time per recipient in Ind',sum(Ind_wait_time_per_round)/np.mean(Avg_Just_Ind_sol))

print ('Waiting time for different blood groups in DDIC',[DDIC_BG_wait[i]/DDIC_BG_matched[i] for i in range(4)],'O,A,B,AB')    
print ('Waiting time for different blood groups in Ind',[Ind_BG_wait[i]/Ind_BG_matched[i] for i in range(4)],'O,A,B,AB')    

print ('Effective DO probability in DDIC',np.mean(Effective_DO_prob_DDIC)*0.2)
print ('Effective DO probability in Ind',np.mean(Effective_DO_prob_Ind)*0.2)
print ('Avg DDIC DO per round',np.mean(DDIC_DO_per_round))
print ('Avg Ind DO per round',np.mean(Ind_DO_per_round)) 

xpos=[i for i in range(1,Runs+1)]
plt.plot(xpos,DDIC_sol_per_round,'r--',label='DDIC solution')
plt.plot(xpos,Ind_sol_per_round,'b--',label='Individual solution')
plt.plot(xpos,DDIC_DO_per_round,'r',label='DDIC dropouts')
plt.plot(xpos,Ind_DO_per_round,'b',label='Individual dropouts')
plt.xlabel('Number of Rounds \n P-U(1,5), PWL-U(1,5), DD-U(1,5), DP=0.2')
plt.ylabel('Number of recipients')
#plt.title('Comparison of Individual solution vs Global solution for 2 countries')
#plt.figure(facecolor="white")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=1.)
plt.savefig('DDIC_vs_Ind.png')