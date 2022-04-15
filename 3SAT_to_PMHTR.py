import sys


def read_file(filename):
    with open(filename, 'r') as f:
        nvars, nclauses = (int(i) for i in f.readline().split())
        sat = []
        for i in range(nclauses):
            sat.append( tuple(int(i) for i in f.readline().split()) )
        return nvars, sat


def print_x_florets(nvars, treefile, lfile):
    for i in range(1, nvars + 1):
        for j in range(1, nvars + 1):
            if i != j:
                treefile.write(f'RR x1_{i}_{j}\n')
                treefile.write(f'x1_{i}_{j} x1l_{i}_{j}\n')
                lfile.write(f'x1l_{i}_{j} R\n')
                treefile.write(f'x1_{i}_{j} x2_{i}_{j}\n')
                treefile.write(f'x2_{i}_{j} x2ll_{i}_{j}\n')
                treefile.write(f'x2_{i}_{j} x2lr_{i}_{j}\n')
                lfile.write(f'x2ll_{i}_{j} {i}\n')
                lfile.write(f'x2lr_{i}_{j} {-i}\n')
                treefile.write(f'x2_{i}_{j} x3_{i}_{j}\n')
                treefile.write(f'x3_{i}_{j} x3ll_{i}_{j}\n')
                treefile.write(f'x3_{i}_{j} x3lr_{i}_{j}\n')
                lfile.write(f'x3ll_{i}_{j} {i}\n')
                lfile.write(f'x3lr_{i}_{j} {-i}\n')
                treefile.write(f'x3_{i}_{j} x4_{i}_{j}\n')
                treefile.write(f'x4_{i}_{j} x4ll_{i}_{j}\n')
                treefile.write(f'x4_{i}_{j} x4lr_{i}_{j}\n')
                lfile.write(f'x4ll_{i}_{j} {j}\n')
                lfile.write(f'x4lr_{i}_{j} {-j}\n')


def print_clause_florets(clauses, treefile, lfile):
    for i in range(len(clauses)):
        treefile.write(f'RR c1_{i}\n')
        treefile.write(f'c1_{i} c5_{i}\n')
        treefile.write(f'c5_{i} c5l_{i}\n')
        lfile.write(f'c5l_{i} C{i}\n')
        # treefile.write(f'c5_{i} c6_{i}\n')
        # treefile.write(f'c6_{i} c6l_{i}\n')
        # lfile.write(f'c6l_{i} C{i}\n')
        treefile.write(f'c1_{i} c2_{i}\n')
        treefile.write(f'c2_{i} c2l_{i}\n')
        lfile.write(f'c2l_{i} C{i}\n')
        treefile.write(f'c2_{i} c3_{i}\n')
        for j in range(3):
            treefile.write(f'c3_{i} c3{j}_{i}\n')
            treefile.write(f'c3{j}_{i} c3{j}ll_{i}\n')
            treefile.write(f'c3{j}_{i} c3{j}lr_{i}\n')
            lfile.write(f'c3{j}ll_{i} {clauses[i][j]}\n')
            lfile.write(f'c3{j}lr_{i} {-clauses[i][j]}\n')
            treefile.write(f'c5_{i} c6{j}_{i}\n')
            treefile.write(f'c6{j}_{i} c6{j}ll_{i}\n')
            treefile.write(f'c6{j}_{i} c6{j}lr_{i}\n')
            lfile.write(f'c6{j}ll_{i} {clauses[i][j]}\n')
            lfile.write(f'c6{j}lr_{i} C{i}\n')


if __name__ == '__main__':
    nvars, clauses = read_file(sys.argv[1])
    with open(sys.argv[1]+'.tree', 'w+') as treefile, open(sys.argv[1]+'.labeling', 'w+') as lfile :
        treefile.write(f'RR RRl\n')
        lfile.write(f'RRl R\n')

        print_x_florets(nvars, treefile, lfile)
        print_clause_florets(clauses, treefile, lfile)
    

