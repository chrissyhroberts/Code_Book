import pandas as pd
import sys
import pyarrow
import time

start = time.time()
p = sys.argv[2]
print (p)
df = pd.read_csv(sys.argv[1])
df.to_parquet(p)
end = time.time()
print(end - start)