
for i in `seq 1 100`;
do
    blasr test_sim1.fa test_sim2.fa -maxLCPLength 15 -m 4 > /dev/null 

done
