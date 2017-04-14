import matplotlib.pyplot as plt

leja_sequence_file = open( "leja_sequence.txt", "r")
lines = leja_sequence_file.readlines()
leja_sequence_file.close()

# read the number of Leja points
nb_leja_points = int(lines[0])
print("Sequence with", nb_leja_points, "Leja points:")

# read and store Leja points values
leja_sequence_x = []
leja_sequence_y = []
values = []

for i in range(0,nb_leja_points):
    values.append(float(lines[i+1]))

for v1 in values:
    for v2 in values:
        leja_sequence_x.append(v1)
        leja_sequence_y.append(v2)

#for v1 in values:
#    for v2 in values:
#        print(v1, " ", v2)

# plot the 2d Leja sequence
plt.scatter(leja_sequence_x,leja_sequence_y,s=1)
plt.title('Nuage de points de Leja')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('lejaSequence_2dDistribution.png')
plt.show()
