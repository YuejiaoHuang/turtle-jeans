# Remove 162482at32523 gene from final data set, as the calculation of dN and dS does not work

# delete gene from final gene list
sed -i '' -e "/162482at32523/d" list_final.txt

# delete faulty codeml output files
rm dN_dS_final_data/162482at32523_dS.dist
rm dN_dS_final_data/162482at32523_dN.dist

# delete gene from pure distance matrices output
rm dist_matrices/162482at32523.dist
