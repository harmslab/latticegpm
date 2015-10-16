import svgwrite

ROTATE = {"U":"R", "R":"D", "D":"L", "L":"U"}

class ConfigurationSVG:

    def __init__(self,sequence, configuration, filename="untitled.svg", rotate=0):
        """ Class for drawing an SVG representation of lattice conformation. """
        n = int(rotate/90)
        for i in range(n):
            configuration = "".join([ROTATE[c] for c in configuration])

        self.sequence = sequence
        self.configuration = configuration
        self._filename = filename
        self.drawing = svgwrite.Drawing(self._filename, size=(400,400))
        self.scale = 15
        self.xoffset = 50
        self.yoffset = 50
        self._bonds = {"d":self._dline,
            "u":self._uline,
            "r":self._rline,
            "l":self._lline,
            ".":self._dot,
            " ":None }
        self._draw_sequence(self.sequence, self.configuration)

    @property
    def svg(self):
        """ Return svg"""
        return self.drawing.tostring()

    def save(self):
        """ save svg """
        self.drawing.save()

    @property
    def filename(self):
        return self._filename

    def add(self, stuff):
        """ add to drawing """
        self.drawing.add(stuff)

    def _line(self,coords):
        """ draw line """
        line = self.drawing.line(coords[0], coords[1], stroke=svgwrite.rgb(10,10,16, '%'))
        return line

    def _uline(self,x,y,*args):
        """ x,y are origin coordinates """
        return self.add(self._line(((x,y+5),(x,y-5))))

    def _dline(self,x,y, *args):
        return self.add(self._line(((x,y-5),(x,y+5))))

    def _rline(self, x,y, *args):
        return self.add(self._line(((x-5,y),(x+5,y))))

    def _lline(self,x,y, *args):
        return self.add(self._line(((x+5,y),(x-5,y))))

    def _letter(self,x,y, char):
        xoffset = -6
        yoffset = 5
        letter = self.drawing.text(char, insert=(x+xoffset,y+ yoffset), style="font-size:20px;font-family:Courier")
        return self.add(letter)

    def _dot(self, x,y, *args):
        xoffset = -5
        yoffset = 0
        letter = self.drawing.text('.', insert=(x+xoffset,y+ yoffset), style="font-size:20px;font-family:Courier")
        return self.add(letter)

    def _build_grid(self,sequence, configuration):
        """ """
        moves = {"U":[0,1], "D":[0,-1], "R":[1,0], "L":[-1,0]}

        # find boundaries for drawing
        xmoves, ymoves = [0], [0]

        # Figure out the perimeter of the configuration
        for i in range(len(configuration)):
            move = configuration[i]
            xmoves.append(xmoves[i] + moves[move][0])
            ymoves.append(ymoves[i] + moves[move][1])

        # bounds of configuration
        xmax = max(xmoves)
        xmin = min(xmoves)
        ymax = max(ymoves)
        ymin = min(ymoves)
        xorigin = 0
        yorigin = 0

        # add letters to grid
        xmax2 = xmax + xmax + 3
        xmin2 = xmin + xmin - 2
        ymax2 = ymax + ymax + 3
        ymin2 = ymin + ymin - 2

        # Size of the 2d grid
        xgrid = abs(xmax2 - xmin2)
        ygrid = abs(ymax2 - ymin2)

        # Build an empty grid object
        grid = [[" " for i in range(ygrid)] for j in range(xgrid)]

        # Set dots on grid
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                if i%2==0 and j%2 == 0:
                    grid[i][j] = "."

        # Realign origin to fit in grid
        grid_xorg = xorigin - xmin2
        grid_yorg = yorigin - ymin2

        # initial position
        xpos = grid_xorg
        ypos = grid_yorg

        grid[xpos][ypos] = sequence[0]

        for i in range(len(configuration)):

            let = sequence[i+1]
            bond = configuration[i].lower()
            # Get direction
            step = moves[configuration[i]]

            # Define change in x and y for both edges and letter
            dx = step[0]
            dy = step[1]

            # Add bond
            xpos += dx
            ypos += dy
            grid[xpos][ypos] = bond

            # Add residue
            xpos += dx
            ypos += dy
            grid[xpos][ypos] = let

        return grid

    def _draw_sequence(self, sequence, configuration):
        """ Draw svg of configuration """
        grid = self._build_grid(sequence, configuration)
        # Draw grid in svg
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                # Try if it's a bond
                try:
                    if self._bonds[grid[i][j]] is not None:
                        item = self._bonds[grid[i][j]]
                        item(self.xoffset+self.scale*i, self.yoffset+self.scale*j)

                # Otherwise its a letter
                except KeyError:
                    self._letter(self.xoffset+self.scale*i, self.yoffset+self.scale*j, grid[i][j])
