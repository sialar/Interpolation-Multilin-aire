import matplotlib.pyplot as plt

data_file = open("interpolation_result_1D.txt", "r")
lines = data_file.readlines()
data_file.close()

size = int(lines[0])

data = []
x_tab, y_tab, real_y_tab = [], [], []

for i in range(1,size+1):
	x = float(lines[i].split(" ")[0])
	y = float(lines[i].split(" ")[1])
	real_y = float(lines[i].split(" ")[2])
	data.append((x,y,real_y))

data = sorted(data, key=lambda col: col[0])

for i in range(len(data)):
	x_tab.append(data[i][0])
	y_tab.append(data[i][1])
	real_y_tab.append(data[i][2])

plt.plot(x_tab, y_tab, "b", linewidth=6, label="Interpolation result")
plt.plot(x_tab, real_y_tab, "r", linewidth=2, label="Real function")#, marker="*")
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
