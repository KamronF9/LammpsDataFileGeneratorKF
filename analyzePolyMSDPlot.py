#PM env deepmd3

import pandas as pd
import matplotlib.pyplot as plt

allMSDs = []
for i in range(3):
    df = pd.read_csv(f'{i}MSD.csv')
    allMSDs.append(df)


# plot all MSDs
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
plt.figure(figsize=(6,4))


for i in range(3):
    # print(allMSDs[i])
    # plt.loglog(allMSDs[i]['dt'],allMSDs[i]['msd_c'])
    plt.plot(allMSDs[i]['dt'],allMSDs[i]['msd_c'])

plt.legend(['Pt/O Surface','Bulk','Pt Surface'])  #, loc=2, prop={"size": 20}
# plt.axis('square')
# plt.ylim(1,2e3)
# plt.xlim(1,2e3)
plt.ylabel("MSD ($\\AA^2$)")
plt.xlabel("t (ps)")
# ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')
plt.savefig('allMSDs.pdf')
plt.close()
