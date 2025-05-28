import pandas as pd
import sys
import ntpath

input = sys.argv[1]

output = input.replace("txt","xlsx")

df = pd.read_csv(input, sep='\t', header=0)
df.to_excel(output, index=False)


