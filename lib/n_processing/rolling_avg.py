import sys
import plotly.express as px

data = [0 for i in range(29903)]

if len(sys.argv) < 2:
    sys.stderr.write(
        "usage: python rolling_avg.py <input bedGraph> [comma-separated window sizes]\n"
    )
    sys.exit(1)

with open(sys.argv[1]) as f:
    for line_s in f:
        line = [field.strip() for field in line_s.split("\t")]
        for i in range(int(line[1]) - 1, int(line[2])):
            if data[i] < 0:
                data[i] = 0
            data[i] += int(line[3])

window_sizes = [500, 1000, 2000, 3000, 5000, 7500]  # default window sizes

if len(sys.argv) >= 3:
    window_sizes = [int(val) for val in sys.argv[2].split(",")]

data_avgs = {}

for window_size in window_sizes:

    data_avg = [0 for i in range(len(data))]
    i = 0
    window_len = min(window_size, len(data))
    total = sum(data[0 : window_len - 1])
    for pos in data_avg:
        window_len = min(window_size, len(data) - i)
        if i > 0:
            value_old = data[i - 1]
        else:
            value_old = 0
        value_new = data[i + window_len - 1]
        total -= value_old
        total += value_new
        data_avg[i] = total / float(window_len)
        i += 1

    data_avgs[str(window_size)] = data_avg

fig = px.line(data_avgs)
fig.show()
