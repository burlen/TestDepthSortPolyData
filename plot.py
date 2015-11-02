import numpy as np
import matplotlib.pyplot as plt

#import run_times.py

old_first_point=(2.5862 + 2.55788 + 2.59058)/3.0
old_cell_bound_cen=(3.14134 + 2.99958 + 2.99566)/3.0
old_param_cen=(3.10709 + 3.3351 + 3.12155)/3.0
new_first_point=(0.823937 + 0.820456 + 0.822529)/3.0
new_cell_bound_cen=(1.02036 + 1.02368 + 1.02367)/3.0
new_param_cen=(1.62802 + 1.64475 + 1.6789)/3.0

case=[1,2,3]
case_name=['point','bounds','param']
old_code=np.array([old_first_point, old_cell_bound_cen, old_param_cen])
new_code=np.array([new_first_point, new_cell_bound_cen, new_param_cen])
speed_up=old_code/new_code

print
print ['code'] + list(case_name)
print ['old'] + list(old_code)
print ['new'] + list(new_code)
print ['speed up'] + list(speed_up)
print

fig = plt.figure()
ax = plt.subplot("121")
plt.plot(case, old_code, 'b-', linewidth=2)
plt.plot(case, new_code, 'r-', linewidth=2)
plt.plot(case, old_code, 'bo', markerfacecolor='w', markeredgecolor='b', markeredgewidth=2, markersize=10, label='old')
plt.plot(case, new_code, 'rv', markerfacecolor='w', markeredgecolor='r', markeredgewidth=2, markersize=10, label='new')
ax.set_xticks(case)
ax.set_xticklabels(case_name)#, fontweight='bold')
ax.set_xlim([0.8, 3.2])
plt.xlabel('sort method', fontweight='bold')
plt.ylabel('seconds', fontweight='bold')
plt.title('run time', fontweight='bold')
plt.grid(True,zorder=0)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.946), numpoints=1)
plt.legend(numpoints=1, ncol=2, prop={'weight':'bold', 'size':10})

ax = plt.subplot("122")
plt.bar(case, speed_up, linewidth=2, edgecolor='k', color='y',zorder=3)
plt.xlabel('sort method', fontweight='bold')
plt.ylabel('speed up', fontweight='bold')
plt.title('speed up', fontweight='bold')
ax.set_xticks([1.4, 2.4, 3.4])
ax.set_xticklabels(case_name)#, fontweight='bold')
ax.set_xlim([0.8, 4.0])
plt.grid(True,zorder=0)

plt.subplots_adjust(wspace=0.25)

plt.savefig('depth_sort_speedup.png', dpi=100)
plt.show()
