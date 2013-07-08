
for i in `seq 1 100`;
do
    blasr test_sim5.fa test_sim6.fa -maxLCPLength 15 -m 4 > /dev/null 

done
