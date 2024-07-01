. ./mfc.sh load -c p -m g

acc=64
mode=--gpu
N=1
n=4

#./mfc.sh run ../../AIAA/AIAA$acc/exp/1/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 11_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/exp/2/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 12_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/exp/3/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 13_$acc -j 8 -q embers $mode
./mfc.sh run ../../AIAA/AIAA$acc/exp/4/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 14_$acc -j 8 -q embers $mode -t simulation

#./mfc.sh run ../../AIAA/AIAA$acc/perlin/1/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 21_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/2/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 22_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/3/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 23_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/4/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 - 24_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/41/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 241_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/42/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 242_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/43/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 243_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/44/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 244_$acc -j 8 -q embers $mode

#./mfc.sh run ../../AIAA/AIAA$acc/sin/1/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 31_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/sin/2/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 32_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/sin/3/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 33_$acc -j 8 -q embers $mode
#./mfc.sh run ../../AIAA/AIAA$acc/sin/4/case.py --case-optimization -N $N -n $n -c phoenix -e batch -a gts-sbryngelson3 -w 08:00:00 -# 34_$acc -j 8 -q embers $mode


#./mfc.sh run ../../AIAA/AIAA$acc/exp/1/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/exp/2/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/exp/3/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/exp/4/case.py --case-optimization -N $N -n $n -c phoenix -t post_process

#./mfc.sh run ../../AIAA/AIAA$acc/perlin/1/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/2/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/3/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/4/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/41/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/42/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/43/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/perlin/44/case.py --case-optimization -N $N -n $n -c phoenix -t post_process

#./mfc.sh run ../../AIAA/AIAA$acc/sin/1/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/sin/2/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/sin/3/case.py --case-optimization -N $N -n $n -c phoenix -t post_process
#./mfc.sh run ../../AIAA/AIAA$acc/sin/4/case.py --case-optimization -N $N -n $n -c phoenix -t post_process


