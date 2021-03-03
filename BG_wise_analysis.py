#Blood Group wise analysis of fusion model with bound on chain length as 2.

import math
import random
import numpy as np
import xlrd
from matplotlib import pyplot as plt
from matplotlib import style
from scipy import stats

random.seed(6)

prob_a_a=4.0/211
prob_a_b=40.0/211
prob_a_ab=15.0/211
prob_a_o=1.0/211
prob_b_a=42.0/211
prob_b_b=2.0/211
prob_b_ab=11.0/211
prob_b_o=4.0/211
prob_ab_a=1.0/211
prob_ab_b=1.0/211
prob_ab_ab=1.0/211
prob_ab_o=1.0/211
prob_o_a=33.0/211
prob_o_b=44.0/211
prob_o_ab=10.0/211
prob_o_o=1.0/211

#cummulative probability x
x=[4.0/211,44.0/211,59.0/211,60.0/211,102.0/211,104.0/211,115.0/211,119.0/211,120.0/211,121.0/211,122.0/211,123.0/211,
   156.0/211,200.0/211,210.0/211,1]

#BG_rec=[] #to store generated recipient BG
#BG_don=[] #to store generated donor BG

#Pair=[]

prob_a=23.0
prob_b=33.2
prob_ab=7.7
prob_o=37.1

y=[0.230,0.562,0.639,1]
#DD=[]


avg_pair_matched=[]
avg_DD_chains=[]
set_add_a=[]
set_add_b=[]
set_add_ab=[]
set_add_o=[]

for k in range(50):
    Pair=[]
    BG_rec=[] #to store generated recipient BG
    BG_don=[] #to store generated donor BG
    RecBG_A=0
    RecBG_B=0
    RecBG_AB=0
    RecBG_O=0
    for i in range(200):
        random_num=random.uniform(0,1)
        if random_num>=0 and random_num <=x[0]:
            
            BG_rec.append('A')
            BG_don.append('A')
            Pair.append(("Pair",i))
            RecBG_A=RecBG_A+1
            
        elif random_num>=x[0] and random_num <=x[1]:
                
            BG_rec.append('A')
            BG_don.append('B')
            Pair.append(("Pair",i))
            RecBG_A=RecBG_A+1
            
        elif random_num>=x[1] and random_num <=x[2]:
                

            BG_rec.append('A')
            BG_don.append('AB')
            Pair.append(("Pair",i))
            RecBG_A=RecBG_A+1
            
        elif random_num>=x[2] and random_num <=x[3]:
                
            BG_rec.append('A')
            BG_don.append('O')
            Pair.append(("Pair",i))
            RecBG_A=RecBG_A+1
            
        elif random_num>=x[3] and random_num <=x[4]:
            
            BG_rec.append('B')
            BG_don.append('A')
            Pair.append(("Pair",i))
            RecBG_B=RecBG_B+1
            
        elif random_num>=x[4] and random_num <=x[5]:
            
            BG_rec.append('B')
            BG_don.append('B')
            Pair.append(("Pair",i))
            RecBG_B=RecBG_B+1
        elif random_num>=x[5] and random_num <=x[6]:
            
            BG_rec.append('B')
            BG_don.append('A')
            Pair.append(("Pair",i))
            RecBG_B=RecBG_B+1
        elif random_num>=x[6] and random_num <=x[7]:
            
            BG_rec.append('B')
            BG_don.append('O')
            Pair.append(("Pair",i))
            RecBG_B=RecBG_B+1
        elif random_num>=x[7] and random_num <=x[8]:
            
            BG_rec.append('AB')
            BG_don.append('A')
            Pair.append(("Pair",i))
            RecBG_AB=RecBG_AB+1
        elif random_num>=x[8] and random_num <=x[9]:
            
            BG_rec.append('AB')
            BG_don.append('B')
            Pair.append(("Pair",i))
            RecBG_AB=RecBG_AB+1
        elif random_num>=x[9] and random_num <=x[10]:
          
            BG_rec.append('AB')
            BG_don.append('AB')
            Pair.append(("Pair",i))
            RecBG_AB=RecBG_AB+1
        elif random_num>=x[10] and random_num <=x[11]:
         
            BG_rec.append('AB')
            BG_don.append('O')
            Pair.append(("Pair",i))
            RecBG_AB=RecBG_AB+1
            
        elif random_num>=x[11] and random_num <=x[12]:
            
            BG_rec.append('O')
            BG_don.append('A')
            Pair.append(("Pair",i))
            RecBG_O=RecBG_O+1
            
        elif random_num>=x[12] and random_num <=x[13]:
            
            BG_rec.append('O')
            BG_don.append('B')
            Pair.append(("Pair",i))
            RecBG_O=RecBG_O+1
            
        elif random_num>=x[13] and random_num <=x[14]:
        
            BG_rec.append('O')
            BG_don.append('AB')
            Pair.append(("Pair",i))
            RecBG_O=RecBG_O+1
            
        elif random_num>=x[14] and random_num <=x[15]:
            
            BG_rec.append('O')
            BG_don.append('O')
            Pair.append(("Pair",i))
            RecBG_O=RecBG_O+1
        else:
            print ('No pair is generated')

    count_a=0
    count_b=0
    count_ab=0
    count_o=0
    for i in range(len(Pair)):
        if BG_rec[i]=='A':
            count_a=count_a+1
        if BG_rec[i]=='B':
            count_b=count_b+1
        if BG_rec[i]=='AB':
            count_ab=count_ab+1
        if BG_rec[i]=='O':
            count_o=count_o+1
    
    Matched_pairs=0
    for i in range(len(Pair)):
        for j in range(len(Pair)):
            if (BG_don[i]==BG_rec[j] and BG_don[j]==BG_rec[i]) and ( BG_rec[j]!='Matched' and BG_rec[j]!='Matched'):
                if BG_don[i]==BG_rec[i]:
                    Matched_pairs=Matched_pairs-1
                    #print ('hello')
                BG_don[i]='Matched'
                BG_rec[j]='Matched'
                BG_don[j]='Matched'
                BG_rec[i]='Matched'
                Pair[i]='Matched'
                Pair[j]='Matched'
                #print('Pair',i,'is matched with Pair',j)
                Matched_pairs=Matched_pairs+2
                break

    #for i in Pair:
       # print (i)

    avg_pair_matched.append(Matched_pairs)
    print ('Number of matched pairs',Matched_pairs)

    rem_a=0 #remaining number of recipient with 'A' blood group
    rem_b=0 #remaining number of recipient with 'B' blood group
    rem_ab=0 #remaining number of recipient with 'AB' blood group
    rem_o=0 #remaining number of recipient with 'O' blood group
    
    for i in range(len(Pair)):
        if BG_rec[i]=='A' and BG_rec[i]!='Matched':
            rem_a=rem_a+1
        if BG_rec[i]=='B' and BG_rec[i]!='Matched':
            rem_b=rem_b+1
        if BG_rec[i]=='AB' and BG_rec[i]!='Matched':
            rem_ab=rem_ab+1
        if BG_rec[i]=='O' and BG_rec[i]!='Matched':
            rem_o=rem_o+1
    
    print ('Number of matched recipients with A blood group, total number', (count_a-rem_a), count_a)
    print ('Number of matched recipients with B blood group', (count_b-rem_b) , count_b)
    print ('Number of matched recipients with AB blood group', (count_ab-rem_ab) , count_ab)
    print ('Number of matched recipients with O blood group', (count_o-rem_o) , count_o)

    DD=[]
    count_a_DD=0
    count_b_DD=0
    count_ab_DD=0
    count_o_DD=0
    for j in range(50):
        random_num=random.uniform(0,1)
        #print (random_num),'random num',
        if random_num>=0 and random_num <=y[0]:
            DD.append('A')
            count_a_DD=count_a_DD+1
        elif random_num>=x[0] and random_num <=y[1]:
            DD.append('B')
            count_b_DD=count_b_DD+1
        elif random_num>=x[1] and random_num <=y[2]:
            DD.append('AB')
            count_ab_DD=count_ab_DD+1
        elif random_num>=x[2] and random_num <=y[3]:
            DD.append('O')
            count_o_DD=count_o_DD+1

    #print ("DD",DD)

    count_DD=0

    for i in range(len(DD)):
        for j in range(len(Pair)):
            if DD[i]==BG_rec[j] and BG_rec[j]!='Matched':
                #print ('DD chain exits, DD', i,'is matched with pair', Pair[j])
                count_DD=count_DD+1
                BG_rec[j]='Matched'
                break
        
    avg_DD_chains.append(count_DD)
    print ('Number of additional benefit through DD', count_DD)
    print ('Number of AB recipient',RecBG_AB)

    add_a=0 #additional number of recipient with 'A' blood group
    add_b=0 #additional number of recipient with 'B' blood group
    add_ab=0 #additional number of recipient with 'AB' blood group
    add_o=0 #additional number of recipient with 'O' blood group
    
    for i in range(len(Pair)):
        if BG_rec[i]=='A' and BG_rec[i]!='Matched':
            add_a=add_a+1
        if BG_rec[i]=='B' and BG_rec[i]!='Matched':
            add_b=add_b+1
        if BG_rec[i]=='AB' and BG_rec[i]!='Matched':
            add_ab=add_ab+1
        if BG_rec[i]=='O' and BG_rec[i]!='Matched':
            add_o=add_o+1

    # print ('Number of additional matched recipients with A blood group', (count_a-add_a)- (count_a-rem_a))
    set_add_a.append(((count_a-add_a)- (count_a-rem_a)))
    #print ('Number of additional matched recipients with B blood group', (count_b-add_b)- (count_b-rem_b))
    set_add_b.append(((count_b-add_b)- (count_b-rem_b)))
    #print ('Number of additional matched recipients with AB blood group', (count_ab-add_ab)- (count_ab-rem_ab))
    set_add_ab.append(((count_ab-add_ab)- (count_ab-rem_ab)))
    #print ('Number of additional matched recipients with O blood group', (count_o-add_o)- (count_o-rem_o))
    set_add_o.append(((count_o-add_o)- (count_o-rem_o)))

print ('Avg no of paired matched',np.mean(avg_pair_matched))
print ('Avg no of additional recipient matched, minimum and maximum ',np.mean(avg_DD_chains), min(avg_DD_chains), max(avg_DD_chains))

print ('set of additional recepient matched for bg A',set_add_a)
print ('set of additional recepient matched for bg B',set_add_b)
print ('set of additional recepient matched for bg AB',set_add_ab)
print ('set of additional recepient matched for bg O',set_add_o)

s_sqr=0
print ('DD chains', avg_DD_chains)

for i in range(len(avg_DD_chains)):
    s_sqr=s_sqr+(avg_DD_chains[i]-np.mean(avg_DD_chains))**2

S=np.sqrt(s_sqr/(len(avg_DD_chains)-1))
#t=(np.mean(avg_DD_chains)-15)/np.sqrt((S)/np.sqrt(len(avg_DD_chains)))
#print ('t statistics',t)

CI_Lower=np.mean(avg_DD_chains)-1.96*(S/np.sqrt(len(avg_DD_chains)))
CI_Upper=np.mean(avg_DD_chains)+1.96*(S/np.sqrt(len(avg_DD_chains)))

print ('Confidence Interval for avg gain with 95% confidence',np.mean(avg_DD_chains), CI_Lower,CI_Upper)

s_sqr=0

for i in range(len(set_add_a)):
    s_sqr=s_sqr+(set_add_a[i]-np.mean(set_add_a))**2

S=np.sqrt(s_sqr/(len(set_add_a)-1))
CI_Lower=np.mean(set_add_a)-1.96*(S/np.sqrt(len(set_add_a)))
CI_Upper=np.mean(set_add_a)+1.96*(S/np.sqrt(len(set_add_a)))

print ('Confidence Interval for A Blood group with 95% confidence',np.mean(set_add_a), CI_Lower,CI_Upper)

s_sqr=0

for i in range(len(set_add_b)):
    s_sqr=s_sqr+(set_add_b[i]-np.mean(set_add_b))**2

S=np.sqrt(s_sqr/(len(set_add_b)-1))
CI_Lower=np.mean(set_add_b)-1.96*(S/np.sqrt(len(set_add_b)))
CI_Upper=np.mean(set_add_b)+1.96*(S/np.sqrt(len(set_add_b)))

print ('Confidence Interval for B Blood group with 95% confidence',np.mean(set_add_b), CI_Lower,CI_Upper)

s_sqr=0

for i in range(len(set_add_ab)):
    s_sqr=s_sqr+(set_add_ab[i]-np.mean(set_add_ab))**2

S=np.sqrt(s_sqr/(len(set_add_ab)-1))
CI_Lower=np.mean(set_add_ab)-1.96*(S/np.sqrt(len(set_add_ab)))
CI_Upper=np.mean(set_add_ab)+1.96*(S/np.sqrt(len(set_add_ab)))

print ('Confidence Interval for AB Blood group with 95% confidence',np.mean(set_add_ab), CI_Lower,CI_Upper)

s_sqr=0

for i in range(len(set_add_o)):
    s_sqr=s_sqr+(set_add_o[i]-np.mean(set_add_o))**2

S=np.sqrt(s_sqr/(len(set_add_o)-1))
CI_Lower=np.mean(set_add_o)-1.96*(S/np.sqrt(len(set_add_o)))
CI_Upper=np.mean(set_add_o)+1.96*(S/np.sqrt(len(set_add_o)))

print ('Confidence Interval for O Blood group with 95% confidence',np.mean(set_add_o), CI_Lower,CI_Upper)


x=[]
for i in range(1,len(avg_pair_matched)+1):
    x.append(i)
y=avg_DD_chains

plt.plot(x,y, linestyle='-')
plt.title('Relative Gain of Fusion model over Independent PKE and DD allocation for 50 replications')
plt.ylabel('Number of additional PKE recipient matched')
plt.xlabel('Number of replications')
plt.show()










           
