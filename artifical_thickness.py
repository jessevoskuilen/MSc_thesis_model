def compute_average_retention(num_slabs):
    total_retentions = [0]  # Total retention after each slab
    total_passing = [0]
    slab_retentions = [0]
    retention_rate = 0.416

    for i in range(1, num_slabs + 1):
        total = retention_rate+retention_rate*total_passing[-1]
        total_passing.append((1-retention_rate)+(1-retention_rate)*total_passing[-1])
        total_retentions.append(total)

    average_retention = sum(total_retentions) / num_slabs
    return [r * 100 for r in slab_retentions], average_retention * 100

# Display values for up to 20 slabs
for n in range(1, 21):
    slab_rets, avg_ret = compute_average_retention(n)
    print(f"Slabs: {n:2d} | Average Retention per slab: {avg_ret:.4f}%")
