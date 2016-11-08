subset2['R23'] = (subset2['O2-3727'] + subset2['O3-4959'] + subset2['O3-5007']) / subset2['Hb']
subset2.plot(x='O3-4363', y='R23', kind='scatter')
plt.show()
