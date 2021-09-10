import sys

bed1 = open(sys.argv[1],"r")
bed2 = open(sys.argv[2],"r")

lines1=bed1.readlines()
lines2=bed2.readlines()


bed1.close()
bed2.close()

#print (len(lines1), len(lines2))

for i in range(0,len(lines1)):
    line1=lines1[i]
    line1 = line1.split()
    chro=line1[0]
    start=line1[1]
    end=line1[2]
    reads1=int(line1[3])
    line2 = lines2[i]
    line2 = line2.split()
    reads2=int(line2[3])
    print(chro + "\t" + start + "\t" + end + "\t" + str(reads1+reads2))
