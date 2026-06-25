import numpy as np
import matplotlib.pyplot as plt

freq_file = "Al2O3.freq"

q_dist = []
freqs = []

with open(freq_file, "r") as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip().split()

    # QE matdyn freq format often has:
    # qx qy qz distance
    # then frequencies over following lines
    if len(line) == 4:
        try:
            qx, qy, qz, x = map(float, line)
            vals = []
            i += 1

            while i < len(lines):
                parts = lines[i].strip().split()

                # next q-point starts
                if len(parts) == 4:
                    try:
                        float(parts[0]); float(parts[1]); float(parts[2]); float(parts[3])
                        break
                    except ValueError:
                        pass

                for p in parts:
                    try:
                        vals.append(float(p))
                    except ValueError:
                        pass

                if len(vals) >= 30:
                    break

                i += 1

            if len(vals) >= 30:
                q_dist.append(x)
                freqs.append(vals[:30])
        except ValueError:
            i += 1
    else:
        i += 1

q_dist = np.array(q_dist)
freqs = np.array(freqs)

print("Number of q-points parsed:", len(q_dist))
print("Number of modes:", freqs.shape[1] if freqs.size else 0)
print("Minimum frequency:", np.min(freqs) if freqs.size else None)

plt.figure(figsize=(8, 5))

for m in range(freqs.shape[1]):
    plt.plot(q_dist, freqs[:, m], linewidth=1)

# high-symmetry vertical lines
labels = [r"$\Gamma$", "X", "K", r"$\Gamma$", "Z", "U", "H", "Z"]
nseg = 7
points_per_segment = 30

tick_indices = [0]
for s in range(1, nseg + 1):
    tick_indices.append(min(s * points_per_segment, len(q_dist)-1))

tick_positions = [q_dist[idx] for idx in tick_indices]

for xpos in tick_positions:
    plt.axvline(xpos, linewidth=0.8)

plt.xticks(tick_positions, labels)
plt.ylabel(r"Frequency (cm$^{-1}$)")
plt.xlabel("Wave vector")
plt.title(r"Al$_2$O$_3$ phonon dispersion")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("Al2O3_phonon_dispersion.png", dpi=300)
plt.show()
