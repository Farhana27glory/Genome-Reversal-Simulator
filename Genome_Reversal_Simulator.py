#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 21:12:27 2020

@author: gloryfz
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########
###########  Random genome reversal simulator for transforming one genome 
###########  into another
###########
###########  Farhana Zaman Glory
###########  with intergenic regions Version 3 with one output file for result dumping,
###########  gene numbers in each reversal and also generating input with 0's as
###########  intergenic genes, this version also tests iterations for each 10 thousand 
###########  with the start and ending time

"""
Created on Thu Feb  6 17:33:57 2020

@author: gloryfz
"""
import math
import re
import random
import sys
import time


num_of_paths = 200
length_limit = 10000
SEEKING_LIMIT = 40

# set the probability for taking three types of reversals (good, neutral, or bad, add up to 1.0)
PROB_BAD, PROB_GOOD = 0.0001, 0.99
#PROB_GOOD = 0.99
prob_bad = PROB_BAD



def intergenic_conversion(Lis):
  for k in range(math.floor(len(Lis))):
    for i in range(len(Lis)):
       if type(Lis[i]) == str:
          M = int(Lis[i])
          Lis.pop(i)
          if M!=0:
            for j in range(M):
              Lis.insert(j+i,0)
          else:
              break
       else:
          continue

  return(Lis)     



#### converting list from string to interger###########

def str_list_to_int_list(str_list):
    int_list = [int(n) for n in str_list]
    return int_list


def Input_Processing(file1,file2):

  f1 = open(file1, 'r') 
  InputA = f1.read()     
  f1.close()              
   ##print (InputA) 0.8


  f2 = open(file2, 'r') 
  InputB = f2.read()     
  f2.close()              
   ##print (InputB) 


  ListA = list(InputA.split(" "))
  ListB = list(InputB.split(" "))
###print(ListA)
###print(ListB)

  ListA =[re.sub('[^-^a-zA-Z0-9]+', '', _) for _ in ListA]
  ListB =[re.sub('[^-^a-zA-Z0-9]+', '', _) for _ in ListB]


  
# using naive method to 
# perform conversion 
  
  
  ListA = [int(ListA[i]) for i in range(0, len(ListA)) if (i % 2==0)]
  ListB = [int(ListB[i]) for i in range(0, len(ListB)) if (i % 2==0)]
  

  ListA =(intergenic_conversion(ListA))
  ListB =(intergenic_conversion(ListB)) 

#### converting list from string to interger###########

  ListA = (str_list_to_int_list(ListA))
  ListB = (str_list_to_int_list(ListB))
  return(ListA, ListB)


ListA, ListB = Input_Processing('testlong4b1.txt','testlong4a.txt')
#ListN = ListB



########### reversal function with sign #############
#@numba.jit(nopython=True)
def reverse_sign(l, start, end): 
    for i in range(start, end): 
         l[i] *= (-1)
         i=i+1
    return l



'''
def reverse(l, first=0, last=-1):
 l = np.array(l) 
 if (first >= last): return
 #reverse_sign(l, first, last)   
 l[first], l[last] = l[last], l[first]
 reverse(l, first+1, last-1)
'''

def reversal_function(l,s,r):
  first_index = 0
  last_index = len(l)-1
  if (s < r):
    if ( first_index < s < last_index and r!=last_index):
        first_half = l[:s]
        #print(first_half)

        middle_half = l[s:r+1]
       #print(middle_half)

        last_half = l[r+1:]
        #print(last_half)

        reverse_sign(middle_half, 0, len(middle_half))
        middle_half.reverse()

        l = first_half + middle_half + last_half
        
    
    elif( s > first_index and r == last_index):
        first_half = l[:s]
        #print(first_half)


        last_half = l[s:r+1]
        #print(last_half)

        reverse_sign(last_half, 0, len(last_half))
        last_half.reverse()

        l = first_half + last_half
        
    
    elif (s == first_index and r != last_index):
        first_half = l[s:r+1]
        #print(first_half)


        last_half = l[r+1:]
        #print(last_half)

        reverse_sign(first_half, 0, len(first_half))
        first_half.reverse()

        l = first_half + last_half
        

    else:
        reverse_sign(l, 0, len(l))
        l.reverse()
        
        
  elif (s == r):
    if(first_index < s < last_index):
        first_half = l[:s]
        #print(first_half)

        middle_half = l[s:r+1]
        #print(middle_half)

        last_half = l[r+1:]
        #print(last_half)

        reverse_sign(middle_half, 0, len(middle_half))
        middle_half.reverse()

        l = first_half + middle_half + last_half
        
        
    elif (s == r == first_index):
        first_half = l[s:r+1]
        #print(first_half)


        last_half = l[r+1:]
        #print(last_half)
    
        reverse_sign(first_half, 0, len(first_half))
        first_half.reverse()

        l = first_half + last_half
        
        
    else:
        first_half = l[:s]
        #print(first_half)


        last_half = l[s:r+1]
        #print(last_half)

        reverse_sign(last_half, 0, len(last_half))
        last_half.reverse()

        l = first_half + last_half
        
  else:
    s, r = r, s
    
    
    if ( first_index < s < last_index and r!=last_index):
        first_half = l[:s]
        #print(first_half)

        middle_half = l[s:r+1]
        #print(middle_half)

        last_half = l[r+1:]
        #print(last_half)

        reverse_sign(middle_half, 0, len(middle_half))
        middle_half.reverse()

        l = first_half + middle_half + last_half
        
        
    elif( s > first_index and r == last_index):
        first_half = l[:s]
        #print(first_half)


        last_half = l[s:r+1]
        #print(last_half)

        reverse_sign(last_half, 0, len(last_half))
        last_half.reverse()

        l = first_half + last_half
        
    
    elif (s == first_index and r != last_index):
        first_half = l[s:r+1]
        #print(first_half)


        last_half = l[r+1:]
        #print(last_half)

        reverse_sign(first_half, 0, len(first_half))
        first_half.reverse()

        l = first_half + last_half
        

    else:
        reverse_sign(l, 0, len(l))
        l.reverse()
        
  return (l)


##### creating another list containing only zero#####
 
Zero = [0 for i in range(len(ListB) + 1)]
Zero = (str_list_to_int_list(Zero))
###print(Zero)




######## creating reversed list with negative signs #########

ListC = [x1 - x2 for (x1, x2) in zip(Zero, ListB)]
###print(ListC)
ListC.reverse()
###print(ListC)




##### adjacency list creation of genome B and genome A ######

def Adjacency_Tracking_Genome(ListB):
 
 ListB = [num for num in ListB if num]
 ListB.insert(len(ListB),ListB[0])
 C = [ListB[i:i+2] for i in range(len(ListB) - 1)]
 ListB = ListB[:-1]
 return C 



# the method for determining the position of the breakpoints
#@numba.jit(nopython=True)
def getBreakPoints(ListA,ListB,ListC):
  listA = ListA
  listB = ListB
  listC = ListC
  breakPoints = [listA[s] for s in range(len(listA)) if ((not (listA[s] in listB)) and (not (listA[s] in listC)))]
  return breakPoints


###### common adjacencies in both genomes #############
#@numba.jit(nopython=True)
def Adjacency_Matching(a,b):
 D = [k for k in a for v in b if (k == v)]
 #D = map(k, a)          
 return len(D)



#@numba.jit(nopython=True)
def Reversal_Simulator():
    
    ListA, ListB = Input_Processing('testlong4b1.txt','testlong4a.txt')
    #ListN = ListB
    ListM = ListA
    ########## ListN now stores reversed ListB as backup ##############
    #reverse_sign(ListN,0,(len(ListN)-1)) 
    #ListN = reversal_function(ListN,0,(len(ListN)-1))
    ###print(ListN)
    
    ######### random reversal generator process ##############
    
    #ListM, ListB = Input_Processing('testaaa.txt','testbbb.txt')
    path_num = 0
    # tries = 0
    # failure = 0
    
    result_file = open('RESULT1.txt', 'w')

    # while (path_num<num_of_paths and tries < tries_limit):
    while (path_num<num_of_paths):
      #ListA = ListM.copy()
      ListA = list(ListM)
      adjacency_listA = Adjacency_Tracking_Genome(ListA)
      adjacency_listB = Adjacency_Tracking_Genome(ListB)
      adjacency_listC = Adjacency_Tracking_Genome(ListC)

      prev_adjacency_count_AB = Adjacency_Matching(adjacency_listA,adjacency_listB)
      prev_adjacency_count_AC = Adjacency_Matching(adjacency_listA,adjacency_listC)
      prev_adjacency_count = max(prev_adjacency_count_AB,prev_adjacency_count_AC)

      reversal_count = 0
      # tries = 0
      # good_taken = 0
      # neutral_taken = 0
      # bad_taken = 0

      while (reversal_count < length_limit):
        # randomly generate a reversal
        # a dice random in range [0.0,1.0)
        # prob_good = PROB_GOOD
        # prob_neutral = 1 - prob_good
        # split_neutral_good = prob_neutral
        # =========================================== pre-assign the prob ===========================================
        dice = random.random()

        adjacency_listA = Adjacency_Tracking_Genome(ListA)
        tmp_reversal_count = reversal_count
        breakPoints = getBreakPoints(adjacency_listA,adjacency_listB,adjacency_listC)
        prob_neutral = (1-PROB_GOOD-PROB_BAD)*((len(adjacency_listB)-len(breakPoints))/len(adjacency_listB))
        # prob_good = 1 - prob_neutral - PROB_BAD
        split_neutral_good = PROB_BAD + prob_neutral
        split_bad_neutral = PROB_BAD
        # split_neutral_good = prob_neutral
        # print(prob_good)
        # dynamic seeking limit (escape function)
        seek_failures = 0
        seeking_limit = SEEKING_LIMIT*(len(breakPoints)/len(adjacency_listB))+20

        while(reversal_count == tmp_reversal_count):
          # print(ListA,dice)
          # print(len(adjacency_listB))
          if (dice >= split_neutral_good):
            # print(breakPoints,len(breakPoints))
            s_breakPt_pos = random.randint(0,len(breakPoints)-1)
            r_breakPt_pos = random.randint(0,len(breakPoints)-1)
            # print(ListA)
            if (s_breakPt_pos == r_breakPt_pos):
              s = ListA.index(breakPoints[s_breakPt_pos][random.choice([0,1])])
              r = ListA.index(breakPoints[s_breakPt_pos][random.choice([0,1])])
            else:
              if (s_breakPt_pos > r_breakPt_pos):
                #tmp_breakPt_pos = s_breakPt_pos
                s_breakPt_pos, r_breakPt_pos = r_breakPt_pos, s_breakPt_pos
                #r_breakPt_pos = tmp_breakPt_pos
              s = ListA.index(breakPoints[s_breakPt_pos][1])
              r = ListA.index(breakPoints[r_breakPt_pos][0])

            #reverse_sign(ListA,s,r)
            ListA = reversal_function(ListA,s,r)
            #print(ListA)
            adjacency_listA = Adjacency_Tracking_Genome(ListA)
            adjacency_count_AB = Adjacency_Matching(adjacency_listA,adjacency_listB)
            adjacency_count_AC = Adjacency_Matching(adjacency_listA,adjacency_listC)
            curr_adjacency_count = adjacency_count_AB + adjacency_count_AC

            if (curr_adjacency_count > prev_adjacency_count):
              prev_adjacency_count = curr_adjacency_count
              reversal_count += 1
              # good_taken += 1
              seek_failures = 0
              # print('breakPoints:',breakPoints,'\n')
            else:
              #reverse_sign(ListA,s,r)
              ListA = reversal_function(ListA,s,r)
              #print(ListA)
              seek_failures += 1
              # print('bk1:',s_breakPt_pos,'bk2',r_breakPt_pos,'\n')
              if (seek_failures > seeking_limit-1): 
                # print("switch to neutral")
                dice = split_bad_neutral
                seek_failures = 0

          elif ((dice < split_neutral_good) and (dice >= split_bad_neutral)):
         
            s_breakPt_pos = random.randint(0,len(breakPoints)-1)
            r_breakPt_pos = random.randint(0,len(breakPoints)-1)
            # print(ListA)
            while(s_breakPt_pos == r_breakPt_pos):
              s_breakPt_pos = random.randint(0,len(breakPoints)-1)
              r_breakPt_pos = random.randint(0,len(breakPoints)-1)
              # s = ListA.index(breakPoints[s_breakPt_pos][0])
              # r = ListA.index(breakPoints[s_breakPt_pos][1])

            if (s_breakPt_pos > r_breakPt_pos):
              #tmp_breakPt_pos = s_breakPt_pos
              #s_breakPt_pos = r_breakPt_pos
              #r_breakPt_pos = tmp_breakPt_pos
              s_breakPt_pos, r_breakPt_pos = r_breakPt_pos, s_breakPt_pos
            s = ListA.index(breakPoints[s_breakPt_pos][1])
            r = ListA.index(breakPoints[r_breakPt_pos][0])

            #reverse_sign(ListA,s,r)
            ListA = reversal_function(ListA,s,r)
            #print(ListA)
            # =======================================================================

            adjacency_listA = Adjacency_Tracking_Genome(ListA)
            adjacency_count_AB = Adjacency_Matching(adjacency_listA,adjacency_listB)
            adjacency_count_AC = Adjacency_Matching(adjacency_listA,adjacency_listC)
            curr_adjacency_count = adjacency_count_AB+adjacency_count_AC
            # prev_adjacency_count = max(prev_adjacency_count_AB,prev_adjacency_count_AC)

            if (curr_adjacency_count == prev_adjacency_count):
              reversal_count += 1
              # neutral_taken += 1
              seek_failures = 0
            else:
              #reverse_sign(ListA,s,r)
              ListA = reversal_function(ListA,s,r)
              #print(ListA)
              seek_failures += 1
              # print('nfnfnfnfnfnfnf',seek_failures,'\n')
              if (seek_failures > seeking_limit-1): 
                # print("switch to good")
                dice = split_neutral_good
                seek_failures = 0
          else:
            s = random.randint(0,(len(ListA)-1))
            r = random.randint(0,(len(ListA)-1))
            #reverse_sign(ListA,min(s,r),max(s,r))
            ListA = reversal_function(ListA,min(s,r),max(s,r))
            #print(ListA)
            adjacency_listA = Adjacency_Tracking_Genome(ListA)
            adjacency_count_AB = Adjacency_Matching(adjacency_listA,adjacency_listB)
            adjacency_count_AC = Adjacency_Matching(adjacency_listA,adjacency_listC)
            curr_adjacency_count = adjacency_count_AB+adjacency_count_AC

            if (curr_adjacency_count < prev_adjacency_count):
              reversal_count += 1
              # bad_taken += 1
              seek_failures = 0
            else:
              #reverse_sign(ListA,s,r)
              ListA = reversal_function(ListA,s,r)
              #print(ListA)
              seek_failures += 1
              if (seek_failures > seeking_limit-1): 
                dice = split_bad_neutral
                seek_failures = 0


        # breakpoints of the new A
        adjacency_listA = Adjacency_Tracking_Genome(ListA)
        breakPoints = getBreakPoints(adjacency_listA,adjacency_listB,adjacency_listC)

        if (len(breakPoints) == 0):
          #print(path_num,reversal_count)
          print(reversal_count)
          result_file.write(str(reversal_count)+'\n')
          path_num += 1
          break

        # if (reversal_count == length_limit):
        #   failure += 1

        # tries += 1
      
    result_file.close()
    



def main():
    
    start = time.perf_counter()
    #start = time.clock()
    #print(ListA)
    #print(ListB)
    #result_file = open('RESULT.txt', 'w')
    Reversal_Simulator()
    #result_file.close()
    end = time.perf_counter()
    #end = time.clock()
    #print('Paths:',path_num,'Failures:',failure,'Failure rate:',float(failure/(path_num+failure)))
    original = sys.stdout
    #print(end)
    original = sys.stdout
    print("Time elapsed during the calculation in seconds:", end - start)

if __name__ == "__main__":
    main()
