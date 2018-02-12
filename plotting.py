import matplotlib.pyplot as plt
import numpy as np

fname = "../awake-baseline/fi04000m.swp"

data = np.loadtxt(fname)

len = data.shape[0]

print data.shape

# zeroth_entry = []
# first_entry = []
# second_entry = []
# third_entry = []
# fourth_entry = []
# fifth_entry = []
# sixth_entry = []
# seventh_entry = []
# eighth_entry = []
# ninth_entry = []
# tenth_entry = []
# eleventh_entry = []
# twelfth_entry = []
#
# for i in range(0, len):
#     zeroth_entry.append(data[i][0])
#     first_entry.append(data[i][1])
#     second_entry.append(data[i][2])
#     third_entry.append(data[i][3])
#     fourth_entry.append(data[i][4])
#     fifth_entry.append(data[i][5])
#     sixth_entry.append(data[i][6])
#     seventh_entry.append(data[i][7])
#     eighth_entry.append(data[i][8])
#     ninth_entry.append(data[i][9])
#     tenth_entry.append(data[i][10])
#     eleventh_entry.append(data[i][11])
#     twelfth_entry.append(data[i][12])
#
# plt.plot(zeroth_entry, second_entry)
# plt.show()


