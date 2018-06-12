import os
import shutil
import re
seq = list()

file_folder = "/N/u/tlicknac/Carbonate/Upstream_Test/"   						#Location of upstream sequences

seq = os.listdir(file_folder)

print(seq)

for i in seq:
	real_i = file_folder + i
        filepath = "/N/u/tlicknac/Carbonate/MEME_Test/" + i                                		#Location of output
        try:
            	comm1 = "mkdir " + "/N/u/tlicknac/Carbonate/MEME_Test/" + i                          	
        except:
               	pass
        os.system(comm1)
        comm2 = "/N/u/tlicknac/Carbonate/software/meme/install/bin/meme " + real_i + " -dna -oc " + filepath + " -nostatus -time 18000 -maxsize 60000 -nmotifs 5 -minw 5 -maxw 18"   #run my meme, not Carbonate's
        os.system(comm2)
