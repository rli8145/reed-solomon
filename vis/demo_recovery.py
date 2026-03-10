#!/usr/bin/env python3
"""Bar chart showing codeword values, corrupted positions in red, recovered values overlaid."""

import json
import matplotlib.pyplot as plt
import numpy as np

with open("vis/demo_data.json") as f:
    data = json.load(f)

n = data["n"]
k = data["k"]
t = data["t"]
codeword = data["codeword"]
corrupted = data["corrupted"]
error_pos = set(data["error_positions"])

x = np.arange(n)
colors = ["#e74c3c" if i in error_pos else "#3498db" for i in range(n)]

fig, ax = plt.subplots(figsize=(12, 5))

# Corrupted values (what was received)
bars = ax.bar(x, corrupted, color=colors, alpha=0.7, label="received (corrupted)")

# Overlay original codeword as hollow markers at error positions
for i in error_pos:
    ax.plot(i, codeword[i], "D", color="#2ecc71", markersize=10, zorder=5,
            markeredgewidth=2, markerfacecolor="none")

# One marker for the legend
ax.plot([], [], "D", color="#2ecc71", markersize=10, markeredgewidth=2,
        markerfacecolor="none", label="original (correct) value")

ax.set_xlabel("position i")
ax.set_ylabel("value (mod {})".format(data["p"]))
msg_str = data.get("message_str", "")
ax.set_title("Reed-Solomon recovery: \"{}\"  (n={}, k={}, t={} errors corrected)".format(msg_str, n, k, t))
ax.set_xticks(x)
ax.legend()
ax.set_ylim(0, data["p"])

plt.tight_layout()
plt.savefig("vis/demo_recovery.png", dpi=150)
plt.show()
