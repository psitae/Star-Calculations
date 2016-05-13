curve = lambda x: x**3
sum = curve(0) + curve(10)
h = (10) / 1000
for i in range(1, 1000, 2):
    sum += 4 * curve(0 + i * h)
for i in range(2, 999, 2):
    sum += 2 * curve(0 + i * h)
sum = sum * h  / 3
print (sum)

