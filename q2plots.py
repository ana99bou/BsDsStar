import matplotlib.pyplot as plt



mb=1.9257122802734448
md0=0.73483032
md1=0.7458868408203149
md2=0.7567413330078149
md3=0.7673809814453149
md4=0.7774560546875026
md5=0.787656860351565

r=[2.785**2*(mb**2+md0**2-2*mb*md0),2.785**2*(mb**2+md0**2-2*mb*md1),2.785**2*(mb**2+md0**2-2*mb*md2),2.785**2*(mb**2+md0**2-2*mb*md3),2.785**2*(mb**2+md0**2-2*mb*md4),2.785**2*(mb**2+md0**2-2*mb*md5)]

print(r)

x_coords = [r[1], r[2], r[3],r[4],r[5]]
x2_coords=[r[0],r[1],r[2],r[4],r[5]]
a0_coords = [0.3164392089843757, 0.3073040771484382, 0.2998211669921883, 0.2998101806640633,0.28969604492187573]
a1_coords = [0.465849609374999, 0.4591717529296865, 0.45148071289062397,0.43830688476562396,0.4274896240234364]
a2_coords = [0.142341918945312, 0.13689941406249956, 0.13215759277343708,0.13119873046874958,0.1251452636718746]
v_coords = [0.10359436035156251, 0.1011114501953125	, 0.09843872070312501,0.09677001953125,0.09193603515625004]

a0_errors = [0.0037851833493419587, 0.006791161606224229, 0.0039762005799883465,0.004205412767346014,0.005326700510947266]
a1_errors = [0.0059065006267652225, 0.0057527241704552964, 0.005821430620981651,0.011746995107780663,0.009569125614957428]
a2_errors = [0.0017166556854391818, 0.0026351162546076226, 0.0018289481048312972,	0.0019957181072026552,0.002208897628274302]
v_errors = [0.0017511783572122202, 0.0019407292847015043, 0.0022020939011495504,0.0025940767769486203,0.0028204901729288936]

# Create a scatter plot with error bars
plt.errorbar(x_coords, a0_coords, yerr=a0_errors, fmt='o', ecolor='blue', capsize=5, capthick=2, elinewidth=1,label=r'$\widetilde{A}_0$')

plt.errorbar(x2_coords, a1_coords, yerr=a1_errors, fmt='o', ecolor='orange', capsize=5, capthick=2, elinewidth=1,label=r'$\widetilde{A}_1$')

plt.errorbar(x_coords, a2_coords, yerr=a2_errors, fmt='o', ecolor='green', capsize=5, capthick=2, elinewidth=1,label=r'$\widetilde{A}_2$')

plt.errorbar(x_coords, v_coords, yerr=v_errors, fmt='o', ecolor='red', capsize=5, capthick=2, elinewidth=1,label=r'$\widetilde{V}$')

plt.axis((0,12,0,0.5))

# Create a scatter plot
#plt.scatter(x_coords, a0_coords, color='blue')
#plt.scatter(x2_coords, a1_coords, color='red')
#plt.scatter(x_coords, a2_coords, color='green')

# Adding labels to the points
#for (x, y) in zip(x_coords, a0_coords):
#    plt.text(x, y, f'({x}, {y})', fontsize=9, ha='right')

# Adding title and labels
#plt.title('A0')
plt.xlabel(r'$q^2[GeV]^2$',fontsize=15)
plt.ylabel('Lattice Form Factors',fontsize=15)

plt.tick_params(axis='both', which='major', labelsize=14) 
plt.legend()
plt.annotate(r'$\bf{preliminary}$',xy=(0.65,0.03),xycoords='axes fraction',fontsize=15,color='magenta',alpha=1)

# Show plot
#plt.grid(True)
plt.savefig('FF-q^2.pdf',transparent=True,bbox_inches='tight')