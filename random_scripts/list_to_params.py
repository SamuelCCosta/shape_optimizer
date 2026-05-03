data = [   0.53446306,    0.50259374,  243.77396884,  -33.96889496,  261.82802112,
    0.42235164,    0.70059529,  233.99220896,  -43.95489961,  288.09370986,
    0.56719944,    0.34218724,  289.65241533, -106.94255163,  203.10652993,
    0.46513083,    0.8759024,   208.72228596,   71.36566829,  130.41681181]

EllipseBundle_name = 'bundle'
# Loop through the list in steps of 5
for i in range(0, len(data), 5):
    # Get the next 5 items
    chunk = data[i:i+5]
    # Join them with commas and print the formatted string
    params = ", ".join(map(str, chunk))
    print(f"{EllipseBundle_name}.add(Ellipse({params}));")