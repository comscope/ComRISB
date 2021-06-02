import matplotlib.pyplot as plt


'''Simple pdf plot for quick visual check.
'''

colors = ['black', 'red', 'green', 'blue', 'orange','violet',
        'darkred', 'darkgreen', 'navy', 'brown', 'chocolate', 'darkorange',
        'gold', 'olive', 'maroon']

def xy_plot(x_list, y_list, xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    plt.plot(x_list, y_list, "-o")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)
    plt.show()


def xy2_plot(x_list, y_list, pattern_list, label_list,
        xlabel='x', ylabel='y', fsave='test.pdf'):
    plt.figure()
    for x, y, pattern, label in zip(x_list, y_list, pattern_list, label_list):
        plt.plot(x, y, pattern, label=label)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(fsave)
    plt.show()
