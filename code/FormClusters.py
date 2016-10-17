import math
import time
import sys
import pysam

def READ_BAM_STATS(file1):
	
	stats = []	
	fp = open(file1,"r")
	for line in fp:
		stats.append(line)
	fp.close()
	return float(stats[0]), float(stats[1]), float(stats[2])

STAT_FILE = sys.argv[1]
MIN_CLUSTER_SIZE = int(sys.argv[2])
[RDL, MEAN_D, SIG_D] = READ_BAM_STATS(STAT_FILE)
SIG_MULT = int(sys.argv[3])
CLUSTER_D = MEAN_D + SIG_MULT*SIG_D - 2*RDL
BP_MARGIN = int(sys.argv[4])
SIG_BOUND = int(sys.argv[5])

print RDL, MEAN_D, SIG_D, CLUSTER_D   
# IGNORE CLUSTERS THAT HAVE A SMALL NUMBER OF MAPPINGS (IN ANALYSIS)

def CompareCluster(map_type, DiscordCluster, Cl2):

   iPos = Cl2.l_bound
   lPos = Cl2.l_bound
   rPos = Cl2.r_bound
	
   val = 100000 # some large value
   
   if map_type == 0 and (lPos - DiscordCluster.l_bound_O) < CLUSTER_D:
        
        if DiscordCluster.C_type=="01" or DiscordCluster.C_type=="10":
	    val = abs(lPos - DiscordCluster.l_bound_O - (rPos - DiscordCluster.r_bound_O))	
          
        elif DiscordCluster.C_type =="00" or DiscordCluster.C_type == "11":
            val = abs(lPos - DiscordCluster.l_bound_O - (DiscordCluster.r_bound_O - rPos))
        
	if val < 2*SIG_MULT*SIG_D:
            return 1
            
   if map_type == 1 and (iPos - DiscordCluster.l_bound_O) < CLUSTER_D:
        return 1

   return 0

def ChangeBound(Cl1, Cl2, typeC):

     	[l_orient, r_orient, indivOrient] = [int(Cl2.C_type[0]), int(Cl2.C_type[1]), int(Cl2.C_type[0])]
	
	if r_orient != 2:
		map_type = 0
	else:
		map_type = 1

        # Can add different conditions for different cluster categories later. Right now, all fall into type 1.
	if typeC == 1:
			
        	if map_type == 0:
		
			# The first condition is redundant currently, but safe	
			if Cl2.l_bound < Cl1.lmin:
				Cl1.lmin = Cl2.l_bound
			elif Cl2.l_bound > Cl1.lmax:
				Cl1.lmax = Cl2.l_bound
			
			if Cl2.r_bound < Cl1.rmin:
				Cl1.rmin = Cl2.r_bound
			elif Cl2.r_bound > Cl1.rmax:
				Cl1.rmax = Cl2.r_bound

			# $ clean up and combine with above? can't, need track of min and max in all cases.
            		if (not l_orient) and Cl2.l_bound > Cl1.l_bound and Cl2.l_bound!=-1:
                		
                		Cl1.l_bound = Cl2.l_bound
                		                                   
            		elif l_orient and Cl2.l_bound < Cl1.l_bound and Cl2.l_bound!=-1:
                		
                		Cl1.l_bound = Cl2.l_bound
                		
            		if (not r_orient) and Cl2.r_bound > Cl1.r_bound and Cl2.r_bound!=-1:
                	
                		Cl1.r_bound = Cl2.r_bound
                		
            		elif r_orient and Cl2.r_bound < Cl1.r_bound and Cl2.r_bound!=-1:
                		
                		Cl1.r_bound = Cl2.r_bound
                		

        	if map_type == 1:
           
			# The first condition is redundant currently due to sorted comparison, but safe
                        if Cl2.l_bound < Cl1.lmin:
                                Cl1.lmin = Cl2.l_bound
                        elif Cl2.l_bound > Cl1.lmax:
                                Cl1.lmax = Cl2.l_bound
 	
            		if not indivOrient and Cl2.l_bound > Cl1.l_bound and Cl2.l_bound!=-1:
                		Cl1.l_bound = Cl2.l_bound
            		elif indivOrient and Cl2.l_bound < Cl1.l_bound and Cl2.l_bound!=-1:
                		Cl1.l_bound = Cl2.l_bound

        
	return Cl1
            
def calculateMargin(List1):

	for item in List1:
		l_orient = int(item.C_type[0])
		r_orient = int(item.C_type[1])
		
		cl_margin_l = MEAN_D + SIG_BOUND*SIG_D - (item.lmax - item.lmin)
		cl_margin_r = MEAN_D + SIG_BOUND*SIG_D  - (item.rmax - item.rmin)

		cl_margin = min(cl_margin_l, cl_margin_r)
		cl_margin = int(math.ceil(cl_margin))

		if cl_margin <= 0:
			# use small margin
			cl_margin = BP_MARGIN

		ERROR_MARGIN = int(BP_MARGIN/2)
		# for SVSim like test where tolerance is only added to start of SV

		#$test-- remove later. Results v good w/SVSim tally on this.
		#cl_margin = ERROR_MARGIN

		if l_orient == 0:
                        item.l_end = item.l_bound + cl_margin
                        item.l_start = item.l_bound - ERROR_MARGIN

                else:
                        item.l_start = item.l_bound - cl_margin
                        item.l_end = item.l_bound + ERROR_MARGIN

		if r_orient == 0:
			item.r_end = item.r_bound + cl_margin
			item.r_start = item.r_bound - ERROR_MARGIN
		else:
			item.r_start = item.r_bound - cl_margin
			item.r_end = item.r_bound + ERROR_MARGIN

# Class definitions may be local for small classes even though used in different modules
# See ReadDiscordant for variable details
class Cluster(object):

    def __init__(self):
            
            self.l_bound=None
            self.r_bound=None
            self.l_bound_O=None  
            self.r_bound_O=None
            
            self.count=0
            self.C_type=None 
            self.Cmap_type = None 

            self.ltid = -1
            self.rtid = -1

	    # min and max locations by read
            self.lmin = -1
            self.lmax = -1
            self.rmin = -1
            self.rmax = -1
            
	    # final cluster locations	
	    self.l_start = -1 
            self.l_end = -1
            self.r_start = -1
            self.r_end = -1	
			
    def __str__(self):
        return "%s %s %s %s %s %s %s %s %s %s" % (self.count, self.C_type,self.ltid, self.l_start, self.l_end, self.l_end - self.l_start, self.rtid, self.r_start, self.r_end, self.r_end - self.r_start)


if __name__== "__main__":
    
    f1=open("../results/text/All_Discords_P_S.txt","r")
    #f2=open("../results/text/All_Discords_I_S.txt","r")
    f3=open("../results/text/All_Clusters.txt","w")
    f4=open("../results/text/VariantMapInp_P.txt","w")
    f5=open("../results/text/Time_FC_loop.txt","w")
    f6=open("../results/text/Time_FC_line.txt","w")
    
    print "Function start"
    
    DList = []
    counter = 0
    newBlockS = 0
    offset = 0

    for line in f1:

        start_fl = time.clock()

        parsed = line.split()

        temp = Cluster()
        temp.l_bound = int(parsed[2])
        temp.r_bound = int(parsed[4])
	temp.lmin = temp.l_bound
	temp.lmax = temp.lmin + 1
	temp.rmin = temp.r_bound
	temp.rmax = temp.rmin + 1
        temp.C_type = parsed[5]
        temp.ltid = parsed[1]
        temp.rtid = parsed[3]
        currentFrag = parsed[0]
        Claimed = 0

        l_orient = int(temp.C_type[0])
        r_orient = int(temp.C_type[1])

        ld = len(DList)

        start_for = time.clock()
    
        for x in range(ld):
	    
	    CC = CompareCluster(0, DList[ld-x-1], temp)
            if DList[ld-x-1].C_type == temp.C_type and CC and DList[ld-x-1].ltid == temp.ltid and DList[ld-x-1].rtid == temp.rtid: 
		
                DList[ld-x-1].count+=1
                DList[ld-x-1] = ChangeBound(DList[ld-x-1], temp, 1)
                Claimed = 1
                f4.write("%s %s\n" % (currentFrag, offset+ld-x))
                break
                # same fragment can be part of different clusters, but same alignment cannot be
                
        end_for = time.clock()

        f5.write("Time taken for %s for loops: %f\n" %(ld, end_for-start_for))
                    
        if not Claimed:

            temp.l_bound_O = temp.l_bound
            temp.r_bound_O = temp.r_bound
            temp.count = 1
	
	    # downsize list	
            if len(DList) > 0:

	 	if temp.ltid != DList[-1].ltid:
                        newBlockS = 0
			calculateMargin(DList)
			for x in range(len(DList)):
                       		if(DList[x].count > MIN_CLUSTER_SIZE):
                                	f3.write("%s %s\n" %(offset+x+1, DList[x]))
                                	
                    	offset+= len(DList)
                        DList = []
                        DList.append(temp)	

		elif temp.ltid == DList[newBlockS].ltid and (temp.l_bound - DList[newBlockS].l_bound_O) > CLUSTER_D:

                    
		    calculateMargin(DList[:newBlockS])	
                    for x in range(newBlockS):
                        if(DList[x].count > MIN_CLUSTER_SIZE):
				f3.write("%s %s\n" %(offset+x+1, DList[x]))
                   		
                    DList = DList[newBlockS:]
		    offset+= newBlockS
                    newBlockS = len(DList)
                    DList.append(temp)
	
		else:
		    DList.append(temp)	
            else:
                    DList.append(temp)
            
            
            Claimed = 1
            f4.write("%s %s\n" % (currentFrag, offset+len(DList)))

        end_fl = time.clock()
        f6.write("Time taken for one line processing %f\n" %(end_fl-start_fl))

    calculateMargin(DList)
    for x in range(len(DList)):
        if(DList[x].count > MIN_CLUSTER_SIZE):
            f3.write("%s %s\n" %(offset + x + 1, DList[x]))
                
    ld1 = offset + len(DList)          
    newBlockS = 0
    offset=0
                    
    f1.close()
    #f2.close()
    f4.close()

    print ld1, " total clusters."
    print "clusters in total (but not all necessarily significant/listed)."

        
    f3.close()
            
   
       
