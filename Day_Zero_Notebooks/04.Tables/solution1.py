bg_merge['total'] = bg_merge['zody'] + bg_merge['ISM'] + bg_merge['thermal']
bg_merge.plot(x='wave', logy=True)
plt.show()
