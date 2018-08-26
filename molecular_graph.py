import re
from itertools import combinations
from math import sqrt
from operator import itemgetter

atomic_radii = dict(Ac=1.88, Ag=1.59, Al=1.35, Am=1.51, As=1.21, Au=1.50, B=0.83, Ba=1.34, Be=0.35, Bi=1.54, Br=1.21,
                    C=0.68, Ca=0.99, Cd=1.69, Ce=1.83, Cl=0.99, Co=1.33, Cr=1.35, Cs=1.67, Cu=1.52, D=0.23, Dy=1.75,
                    Er=1.73, Eu=1.99, F=0.64, Fe=1.34, Ga=1.22, Gd=1.79, Ge=1.17, H=0.23, Hf=1.57, Hg=1.70, Ho=1.74,
                    I=1.40, In=1.63, Ir=1.32, K=1.33, La=1.87, Li=0.68, Lu=1.72, Mg=1.10, Mn=1.35, Mo=1.47, N=0.68,
                    Na=0.97, Nb=1.48, Nd=1.81, Ni=1.50, Np=1.55, O=0.68, Os=1.37, P=1.05, Pa=1.61, Pb=1.54, Pd=1.50,
                    Pm=1.80, Po=1.68, Pr=1.82, Pt=1.50, Pu=1.53, Ra=1.90, Rb=1.47, Re=1.35, Rh=1.45, Ru=1.40, S=1.02,
                    Sb=1.46, Sc=1.44, Se=1.22, Si=1.20, Sm=1.80, Sn=1.46, Sr=1.12, Ta=1.43, Tb=1.76, Tc=1.35, Te=1.47,
                    Th=1.79, Ti=1.47, Tl=1.55, Tm=1.72, U=1.58, V=1.33, W=1.37, Y=1.78, Yb=1.94, Zn=1.45, Zr=1.56)

atomic_numbers = dict(H='01', He='02', Li='03', Be='04', B='05', C='06', N='07', O='08', F='09', Ne='10', Na='11',
                    Mg='12', Al='13', Si='14', P='15', S='16', Cl='17', Ar='18', K='19', Ca='20', Sc='21', Ti='22',
                    V='23', Cr='24', Mn='25', Fe='26', Co='27', Ni='28', Cu='29', Zn='30', Ga='31', Ge='32', As='33',
                    Se='34', Br='35', Kr='36', Rb='37', Sr='38', Y='39', Zr='40', Nb='41', Mo='42', Tc='43', Ru='44',
                    Rh='45', Pd='46', Ag='47', Cd='48', In='49', Sn='50', Sb='51', Te='52', I='53', Xe='54', Cs='55',
                    Ba='56', La='57', Ce='58', Pr='59', Nd='60', Pm='61', Sm='62', Eu='63', Gd='64', Tb='65', Dy='66',
                    Ho='67', Er='68', Tm='69', Yb='70', Lu='71', Hf='72', Ta='73', W='74', Re='75', Os='76', Ir='77',
                    Pt='78', Au='79', Hg='80', Tl='81', Pb='82', Bi='83', Po='84', At='85', Rn='86', Fr='87', Ra='88',
                    Ac='89', Th='90', Pa='91', U='92', Np='93', Pu='94', Am='95', Cm='96', Bk='97', Cf='98', Es='99',
                    Fm='100', Md='101', No='102', Lr='103', Rf='104', Db='105', Sg='106', Bh='107', Hs='108', Mt='109',
                    Ds='110', Rg='111', Cn='112', Nh='113', Fl='114', Mc='115', Lv='116', Ts='117', Og='118')

class Graph:
    """Represents a molecular graph."""
    __slots__ = ['elements', 'x_coordinates', 'y_coordinates', 'z_coordinates', 'adjacency_list', 'atomic_radii',
                 'atomic_numbers']

    def __init__(self, elements=None, x_coordinates=None, y_coordinates=None, z_coordinates=None,
                 adjacency_list=None, atomic_radii=None, atomic_numbers=None):
        if elements:
            self.elements = elements
        else:
            self.elements = []
        if x_coordinates:
            self.x_coordinates = x_coordinates
        else:
            self.x_coordinates = []
        if y_coordinates:
            self.y_coordinates = y_coordinates
        else:
            self.y_coordinates = []
        if z_coordinates:
            self.z_coordinates = z_coordinates
        else:
            self.z_coordinates = []
        if adjacency_list:
            self.adjacency_list = adjacency_list
        else:
            self.adjacency_list = {}
        if atomic_radii:
            self.atomic_radii = atomic_radii
        else:
            self.atomic_radii = []
        if atomic_numbers:
            self.atomic_numbers = atomic_numbers
        else:
            self.atomic_numbers = []

    def read_file(self, file_path: str) -> None:
        """Reads an XYZ file, searches for elements and their cartesian coordinates
        and adds them to corresponding arrays"""
        pattern = re.compile(r'([A-Za-z]{1,3})\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)')
        with open(file_path) as file:
            for element, x, y, z in pattern.findall(file.read()):
                self.elements.append(element)
                self.x_coordinates.append(float(x))
                self.y_coordinates.append(float(y))
                self.z_coordinates.append(float(z))
        self.atomic_radii = [atomic_radii[element] for element in self.elements]
        self.atomic_numbers = [atomic_numbers[element] for element in self.elements]
        self._generate_adjacency_list()

    def _generate_adjacency_list(self):
        """Generates an adjacency list from atomic cartesian coordinates"""
        node_ids = range(len(self.elements))
        for i, j in combinations(node_ids, 2):
            x_i, y_i, z_i = self.__getitem__(i)[1]
            x_j, y_j, z_j = self.__getitem__(j)[1]
            distance = sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
            if 0.1 < distance < (self.atomic_radii[i] + self.atomic_radii[j]) * 1.3:
                self.adjacency_list.setdefault(i, set()).add(j)
                self.adjacency_list.setdefault(j, set()).add(i)

    def edges(self, adjacency_list=None):
        """Creates an iterator with all graph edges"""
        edges = set()
        for node, neighbours in self.adjacency_list.items():
            for neighbour in neighbours:
                edge = frozenset([node, neighbour])
                if edge in edges:
                    continue
                edges.add(edge)
                yield node, neighbour


    def subgraph(self, nodes, diff=False):
        if not diff:
            adj_list = {node: self.adjacency_list[node] & set(nodes)
                        for node in self.adjacency_list.keys() & set(nodes)}
        else:
            adj_list = {node: self.adjacency_list[node] - set(nodes)
                        for node in self.adjacency_list.keys() - set(nodes)}
        if nodes:
            els = itemgetter(*nodes)(self.elements)
            xs = itemgetter(*nodes)(self.x_coordinates)
            ys = itemgetter(*nodes)(self.y_coordinates)
            zs = itemgetter(*nodes)(self.z_coordinates)
            rs = itemgetter(*nodes)(self.atomic_radii)
            nums = itemgetter(*nodes)(self.atomic_numbers)
        else:
            els = self.elements
            xs = self.x_coordinates
            ys = self.y_coordinates
            zs = self.z_coordinates
            rs = self.atomic_radii
            nums = self.atomic_numbers
        return Graph(elements=els, x_coordinates=xs, y_coordinates=ys, z_coordinates=zs,
                 adjacency_list=adj_list, atomic_radii=rs, atomic_numbers=nums)


    def atom_string(self, node):
        string = self.atomic_numbers[node]
        neighbours = self.adjacency_list[node]
        if neighbours:
            string += ''.join(sorted(itemgetter(*neighbours)(self.atomic_numbers)))
        return string

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, position):
        return self.elements[position], (
            self.x_coordinates[position], self.y_coordinates[position], self.z_coordinates[position])