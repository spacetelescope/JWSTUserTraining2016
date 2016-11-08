df['O3_4363/5007'] = df['O3-4363'] / df['O3-5007']
df.plot(x='O3-5007', y='O3_4363/5007', kind='scatter')
plt.show()
