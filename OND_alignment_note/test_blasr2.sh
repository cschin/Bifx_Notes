
for i in `seq 1 100`;
do
    blasr test_sim3.fa test_sim4.fa -maxLCPLength 15 -m 4 > /dev/null 

done
