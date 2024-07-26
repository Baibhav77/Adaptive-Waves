import pandas as pd

# Read the CSV file, assuming the first row is the header
df = pd.read_csv('plain.csv')

# Write the data back to CSV without the header
df.to_csv('plain_1.csv', index=False, header=False)

print("Headers removed and saved to 'plain_1.csv'.")
