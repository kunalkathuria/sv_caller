# Can include check to see if map type = 1 if other documented end of random insertion was claimed or not. If not, belongs to another variant.
# Use formatBP.py to make bedpe documentation amenable to this reading format

import sys
from collections import Counter

def file_len(f):
    
    for i, l in enumerate(f):
        pass
    return i + 1

if __name__ == "__main__":

	f1 = open("../results/text/All_Variants.txt","r")
        f8 = open("../results/text/DisjSetCover_S.txt","r")
	f11 = open("../results/text/allPositives.txt","w")
	f12 = open("../results/text/allPositives.bedpe","w")
	f13 = open("../results/text/deletions.bedpe","w")
        f14 = open("../results/text/tandemDuplications.bedpe","w")
        f15 = open("../results/text/inversions.bedpe","w")
        f16 = open("../results/text/insertions.bedpe","w")
        f17 = open("../results/text/unknowns.bedpe","w")
		
	DisjSC = []

        for line in f8:
				
                DisjSC.append(int(line))
		
        x = Counter(DisjSC)
	#print x
	#print Claimed
	counter = -1

	for line2 in f1:
		counter+=1
		line2_split = line2.split()
		num = int(line2_split[0])
		if x[num] > 0:
                    
			#print num
			f11.write("%s\n" %line2)
			
			# If not insertion as bp3 is not set
			if line2_split[1] == "DEL":
                            
                            f13.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1]))
                            
                        elif line2_split[1][0:2] == "TD":

			    [bp1_s, bp1_e] = min(int(line2_split[3]),int(line2_split[6])), min(int(line2_split[4]), int(line2_split[7]))
			    [bp2_s, bp2_e] = max(int(line2_split[3]),int(line2_split[6])), max(int(line2_split[4]), int(line2_split[7]))
 
			    f14.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], bp1_s, bp1_e, line2_split[5], bp2_s, bp2_e,line2_split[1]))				

			elif line2_split[1] == "INV":
                            f15.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1]))

                        elif line2_split[1][0:2] == "IN":

			    # two lines for insertion as in bedpe format
			    [bp2_s, bp2_e] = min(int(line2_split[6]),int(line2_split[9])), min(int(line2_split[7]), int(line2_split[10]))
                            [bp3_s, bp3_e] = max(int(line2_split[6]),int(line2_split[9])), max(int(line2_split[7]), int(line2_split[10]))

                            f16.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(line2_split[5], bp2_s, bp3_s, line2_split[2], line2_split[3], line2_split[4], line2_split[1]))

                        # this is read from cluster file, so has 2 bps only
			elif line2_split[1] == "Unknown":

                            f17.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %(line2_split[2], line2_split[3], line2_split[4], line2_split[5], line2_split[6], line2_split[7],line2_split[1]))
                            
	    
	
	f1.close()
	f8.close()
        f11.close()
	f12.close()
	f13.close()
	f14.close()
	f15.close()
	f16.close()
	f17.close()




                                        
                        
