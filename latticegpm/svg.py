import svgwrite

ROTATE = {"U":"R", "R":"D", "D":"L", "L":"U"}

class Add:
    
    def __init__(self, drawing, fontsize=20):
        """ Class for adding objects to lattice model svgwrite drawing object."""
        self.drawing = drawing
        self.fontsize = fontsize
        self.stepsize = 0.25*fontsize
        self.linewidth = 0.1*fontsize
        
    def line(self, start, end):        
        """ Create an svg line object to add to a svgwrite drawing object.
        """
        line = self.drawing.line(start, end, stroke=svgwrite.rgb(10,10,16, '%'), style="stroke-width:" + str(self.linewidth))
        self.drawing.add(line)
    
    def right(self, x, y):
        """ Add a right facing line to figure. """
        start = (x-self.stepsize,y)
        end = (x+self.stepsize,y)
        self.line(start,end)
        
    def left(self, x, y):
        """ Add a right facing line to figure. """
        start = (x+self.stepsize,y)
        end = (x-self.stepsize,y)
        self.line(start,end)

    def up(self, x, y):
        """ Add a right facing line to figure. """
        start = (x,y-self.stepsize)
        end = (x,y+self.stepsize)
        self.line(start, end)

    def down(self, x, y):
        """ Add a right facing line to figure. """
        start = (x,y+self.stepsize)
        end = (x,y-self.stepsize)
        self.line(start,end)

    def letter(self, x, y, char):
        """ Add letter to grid where residues exist. """
        xoffset = -self.stepsize-0.1*self.stepsize
        yoffset = self.stepsize
        letter = self.drawing.text(char, insert=(x+xoffset,y+ yoffset), style="font-size:"+str(self.fontsize)+"px;font-family:Courier")
        return self.drawing.add(letter)
        
    def dot(self, x, y):
        """ Add a dot on grid where no letter exists"""
        xoffset = -self.stepsize
        yoffset = 0
        letter = self.drawing.text('.', insert=(x+xoffset,y+ yoffset), style="font-size:"+str(self.fontsize)+"px;font-family:Courier")
        return self.drawing.add(letter)

class Box(object):
    
    def __init__(self, sequence, configuration):
        """ Build the svg box object. """
        # Create a grid object to get box size and place.
        self.elements = self.build(sequence, configuration)
        self.height = len(self.elements)
        self.width = len(self.elements[0])

        
    @staticmethod
    def build(sequence, configuration):
        """ """
        moves = {"U":[0,-1], "D":[0,1], "R":[1,0], "L":[-1,0]}

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
        grid = [[" " for i in range(xgrid)] for j in range(ygrid)]

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

        grid[ypos][xpos] = sequence[0]

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
            grid[ypos][xpos] = bond

            # Add residue
            xpos += dx
            ypos += dy
            grid[ypos][xpos] = let

        return grid

class Drawing(svgwrite.Drawing):
    
    def __init__(self, box, filename="untitled.svg", element_size=20):
        """ Hijack svgwrite's Drawing object and fill in with specific settings. """
        
        # Set grid
        self.box = box
        
        # Calculate the ysize of each element
        self.yelement_size = element_size
        
        # Use pixel aspect ratio to resize xsize
        self.xelement_size = element_size
        
        self.height = self.yelement_size * self.box.height
        self.width = self.xelement_size * self.box.width
        self.shape = [int(self.width), self.height]
        self.scale = element_size
        self.offset = 0.25 * self.scale
        
        # Construct drawing
        super(Drawing, self).__init__(filename, size=self.shape)
        
        # Adding class
        self.add_to_svg = Add(self, fontsize=element_size)
        
        # Object for how to draw a configuration
        self.svg_mapping = {"d":self.add_to_svg.down,
            "u":self.add_to_svg.up,
            "r":self.add_to_svg.right,
            "l":self.add_to_svg.left,
            ".":self.add_to_svg.dot,
            " ":None
        }
        self.draw()
    
    def draw(self):
        """ Draw this drawing. """
        # Draw grid in svg
        for y in range(self.box.height):
            
            for x in range(self.box.width):
                
                # Try if it's a bond
                try:
                    
                    if self.svg_mapping[self.box.elements[y][x]] is not None:
                        
                        # Get SVG element to add
                        item = self.svg_mapping[self.box.elements[y][x]]
                        
                        # Add element
                        item(2*self.offset+self.scale*x, 2*self.offset+self.scale*y)
                        
                        #item(self.width+self.scale*i, self.height+self.scale*j)

                # Otherwise its a letter
                except KeyError:
                    
                    self.add_to_svg.letter(2*self.offset+self.scale*x, 2*self.offset+self.scale*y, self.box.elements[y][x])
                    
                    #self.add_to_svg.letter(self.width+self.scale*i, self.height+self.scale*j, self.box.elements[i][j])
                    
        


class Configuration:

    def __init__(self, sequence, configuration, filename="untitled.svg", rotation=0, fontsize=20):
        """ Class for drawing an SVG representation of lattice conformation. """
        
        # Rotate configuration if given
        configuration, rotation = self.rotate(configuration, rotation)
        
        # Set sequence and configuration
        self.sequence = sequence
        self.configuration = configuration
        self.filename = filename
        self.rotation = rotation
        self.fontsize = fontsize
        
        # Build SVG grid object
        self.box = Box(self.sequence, self.configuration)
        self.drawing = Drawing(self.box, filename=self.filename, element_size=self.fontsize)
        
    @property
    def string(self):
        """ Return svg as a string. """
        return self.drawing.tostring()

    @property
    def notebook(self):
        """ Display SVG in Jupyther notebook. """
        
        # Import IPython display for notebook
        from IPython.display import SVG as ipython_display
        
        # Display in notebook
        return ipython_display(self.string)
    
    @staticmethod
    def rotate(configuration, rotation):
        # Rotate svg if desired.
        n = int(rotation/90)
        for i in range(n):
            configuration = "".join([ROTATE[c] for c in configuration])
        return configuration, rotation
    
    def save(self):
        """ save svg """
        self.drawing.save()
        
    def draw(self):
        """ Draw svg of configuration """   
        # If using IPython, try to print to notebook
        try:
            return self.notebook
            
        except ImportError:
            raise Warning(""" IPython not installed. """)