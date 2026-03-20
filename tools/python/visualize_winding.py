import matplotlib.pyplot as plt

class Box:
    def __init__(self, scale, vc="red", ec="blue", offset=0):
        self.scale = scale
        self.offset = offset
        self.verticies = [[scale, -scale, -scale], [scale, scale, -scale], [-scale, -scale, -scale], [-scale, scale, -scale], [scale, -scale, scale], [scale, scale, scale], [-scale, -scale, scale], [-scale, scale, scale]]
        self.edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [3, 7], [3, 2], [2, 6], [4, 5], [4, 6], [5,  7], [7, 6]] 
        self.edge_color = ec
        self.vertex_color = vc
         
    def plot(self, ax, vlabel=True, elabel=False):
        for vertID, vert in enumerate(self.verticies):
            ax.scatter(vert[0], vert[1], vert[2], c=self.vertex_color)
            if vlabel:
                ax.text(vert[0], vert[1], vert[2], f"{self.offset + vertID}", fontsize=25)
        for edge in self.edges:
            ax.plot([self.verticies[edge[0]][0], self.verticies[edge[1]][0]], [self.verticies[edge[0]][1], self.verticies[edge[1]][1]], [self.verticies[edge[0]][2], self.verticies[edge[1]][2]], color=self.edge_color)

class Wedge:
    def __init__(self, A, B, ec="green"):
        self.A = A
        self.B = B
        self.edge_color = ec
    def plot(self, ax):
        for vA, vB in zip(self.A.verticies, self.B.verticies):
            ax.plot([vA[0], vB[0]], [vA[1], vB[1]], [vA[2], vB[2]], color=self.edge_color)

def main():
    core = Box(0.5)
    envelope = Box(2, offset=8)
    star = Wedge(core, envelope)
    infinity = Box(5, offset=16)
    vacuum = Wedge(envelope, infinity)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={"projection": "3d"})
    core.plot(ax)
    envelope.plot(ax)
    star.plot(ax)
    infinity.plot(ax)
    vacuum.plot(ax)
    ax.view_init(30, 30)
    plt.show()

if __name__ == "__main__":
    main()
