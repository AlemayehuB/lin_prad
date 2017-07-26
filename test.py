import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

fig  = plt.figure()
fig.set_figwidth(6.0)
fig.set_figheight(6.0)
ax = fig.add_subplot(1,1,1)
ax.pcolormesh([1,2,3,4],[1,2,3,4],[10,12,17,19])
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
ax.set_title("Proton Locations")
plt.show()
