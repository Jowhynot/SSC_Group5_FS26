import pandas as pd
import matplotlib.pyplot as plt

file_path = r"C:\Users\jsjan\Downloads\SSC Data&Codes\G5_mol film\G5\G5_molten_film_temperature.dat"

# Read whitespace-separated data
df = pd.read_csv(
    file_path,
    sep=r"\s+",
    header=None,
    names=["seconds", "temperature"]
)

plt.figure(figsize=(10, 5))
plt.plot(df["seconds"], df["temperature"])

plt.xlabel("Seconds")
plt.ylabel("Temperature")
plt.title("Molten Film Temperature vs Time")
plt.grid(True)
plt.tight_layout()

plt.show()