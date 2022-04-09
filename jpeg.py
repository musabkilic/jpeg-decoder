import binascii
import struct
import time
import sys
from functools import reduce

class Decoder:
    def __init__(self, im):
        self.im = im
        self.j = 0
        self.i = 0
        self.n = len(self.im)

        self.QUANTIZATION_TABLES = {}

    def int2byte(self, x):
        return x.to_bytes(1, sys.byteorder)

    def read(self, x):
        buffer = self.im[self.i:self.i+x]
        self.i += x
        return buffer

    def readb(self, x):
        return [struct.unpack("B", self.int2byte(y))[0] for y in self.read(x)]

    def readh(self, x):
        return [binascii.hexlify(self.int2byte(y)).decode().upper() for y in self.read(x)]

    def log(self, *args, **kwargs):
        text = " ".join(map(str, args))
        print(f"{hex(self.j)[2:]:0>5}", text)
        self.j += kwargs.get("move", 0)

    def pprint_table(self, table):
        for row in table:
            print("["+ ", ".join(["{:3}".format(x) for x in row]) + "]")

    def dezigzag(self, flat):
        table = [[0] * 8 for _ in range(8)]
        x, y = 0, 0
        n, d = 0, 1
        l = 0
        while n < 64:
            table[y][x] = flat[n]
            if l:
                x += -d
                y += +d
                if not ((0 <= x - d <= 7) & (0 <= y + d <= 7)):
                    l = not l
                    d = -d
            else:
                x += (n >= 35) ^ (d == 1)
                y += (n >= 35) ^ (d != 1)
                l = not l
            n += 1
        return table

    def start_of_image(self):
        ### START OF IMAGE
        ### SOI
        self.log(":START OF IMAGE", move=2)

    def start_of_frame(self, n):
        ### START OF FRAME
        ### SOF

        length = struct.unpack(">H", self.read(2))[0]
        self.log(":START OF FRAME", n-0xC0, move=2)

        precision = self.readb(1)[0]
        nlines = struct.unpack(">H", self.read(2))[0]
        nsamples = struct.unpack(">H", self.read(2))[0]
        ncomponents = self.readb(1)[0]

        self.log(":PRECISION", precision, move=1)
        self.log(":NUMBER OF LINES", nlines, move=2)
        self.log(":NUMBER OF SAMPLES PER LINE", nsamples, move=2)
        self.log(":NUMBER OF COMPONENTS", ncomponents, move=1)

        assert ncomponents * 3 + 8 == length

        for i in range(ncomponents):
            cidentifier = self.readb(1)[0]
            hvfactor = self.readb(1)[0]
            hfactor = hvfactor >> 4
            vfactor = hvfactor & 0x0F
            destination = self.readb(1)[0]
            self.log(":PARAMETERS", hfactor, vfactor, destination, move = 3)

        #self.read(length-2)
        #print(" ".join(self.readh(ncomponents*3)))

    def huffman_table(self):
        ### DEFINE HUFFMAN TABLE
        ### DHT

        length = struct.unpack(">H", self.read(2))[0]

        self.log(":HUFFMAN TABLES", length, move=2)
        i = 0
        while i < length - 2:
            bits = []
            table = self.readb(1)[0]
            tclass = table >> 4
            tdestination = table & 0x0F

            for x in range(16):
                bits.append(self.readb(1)[0])

            huffval = self.readb(sum(bits))

            huffsize = []
            I, J = 1, 1
            while I <= 16:
                if J > bits[I - 1]:
                    I += 1
                    J = 1
                else:
                    huffsize.append(I)
                    J += 1

            huffcode = [0]
            for x in range(1, len(huffsize)):
                huffcode.append((huffcode[-1] + 1) << (huffsize[x] - huffsize[x - 1]))

            ehufco = [0] * (max(huffval) + 1)
            ehufsi = [0] * (max(huffval) + 1)
            for x in range(len(huffsize)):
                I = huffval[x]
                ehufco[I] = huffcode[x]
                ehufsi[I] = huffsize[x]

            self.log(":HUFFMAN TABLE", ["DC", "AC"][tclass], tdestination, move=17 + sum(bits))
            #print(bits)
            #print(huffsize)
            print(huffcode)
            #print(ehufco)
            #print(ehufsi)
            #print(" ".join(self.readh(sum(bits))))
            i += 17 + sum(bits)
        assert i == length - 2

        #print(" ".join(self.readh(length - 3)))

    def quantization_table(self):
        ### DEFINE QUANTIZATION TABLE
        ### DQT

        length = struct.unpack(">H", self.read(2))[0]

        self.log(":QUANTIZATION TABLES", length, move=2)
        i = 0
        while i < length - 2:
            table = []
            prec_dest = self.readb(1)[0]
            prec = prec_dest >> 4
            dest = prec_dest & 0x0F

            for x in range(64):
                table.append(struct.unpack("H" if prec else "B", self.read(1))[0])

            self.log(":QUANTIZATION TABLE", prec, dest, move=65 + 64 * prec)
            self.pprint_table(self.dezigzag(table))
            self.QUANTIZATION_TABLES[dest] = table
            i += 65 + 64 * prec
        assert i == length - 2
        self.dezigzag(table)

    def start_of_scan(self):
        ### START OF SCAN
        ### SOS

        length = struct.unpack(">H", self.read(2))[0]
        self.log(":START OF SCAN", length, move=2)

        ncomponents = self.readb(1)[0]
        self.log(":NUMBER OF COMPONENTS", ncomponents, move=1)

        assert length == 6 + 2 * ncomponents

        for i in range(ncomponents):
            cselector = self.readb(1)[0]
            dacentropy = self.readb(1)[0]
            dcentropy = dacentropy >> 4
            acentropy = dacentropy & 0x0F

            self.log(":SCAN COMPONENT", cselector, move=1)
            self.log(":DC - AC ENTROPY", dcentropy, acentropy, move=1)

        sspectral = self.readb(1)[0]
        espectral = self.readb(1)[0]
        hlposition = self.readb(1)[0]
        hposition = hlposition >> 4
        lposition = hlposition & 0x0F

        self.log(":START OF SELECTION", sspectral, move=1)
        self.log(":END OF SELECTION", espectral, move=1)
        self.log(":HIGH - LOW BIT TRANSFORM", hposition, lposition, move=1)

    def application_zero(self):
        ### APPLICATION0
        ### APP0

        length = struct.unpack(">H", self.read(2))[0]
        self.log(":APPLICATION0", move=2)

        identifier = struct.unpack("sssss", self.read(5))
        version = str(self.readb(1)[0]) + "." + str(self.readb(1)[0])
        units = self.readb(1)[0]
        xdensity = struct.unpack(">H", self.read(2))[0]
        ydensity = struct.unpack(">H", self.read(2))[0]
        xthumbnail = self.readb(1)[0]
        ythumbnail = self.readb(1)[0]

        assert length == 16, "THUMBNAILS NOT IMPLEMENTED YET"

        self.log(":LENGTH", length, move=2)
        self.log(f":IDENTIFIER \"{b''.join(identifier).decode()}\"", move=5)
        self.log(":VERSION", version, move=2)
        self.log(":UNITS", units, move=1)
        self.log(":XDENSITY", xdensity, move=2)
        self.log(":YDENSITY", ydensity, move=2)
        self.log(":XTHUMBNAIL", xthumbnail, move=1)
        self.log(":YTHUMBNAIL", ythumbnail, move=1)

    def application_specific(self, n):
        ### APPLICATION SPECIFIC
        ### APPn

        length = struct.unpack(">H", self.read(2))[0]
        self.log(":APPLICATION SPECIFIC", n-0xE0, move=length)

        self.read(length-2)#print("".join(self.readh(length-2)))

    def comment(self):
        ### COMMENT
        ### COM

        length = struct.unpack(">H", self.read(2))[0]
        self.log(":COMMENT", move=2)
        print(b"".join(struct.unpack("s"*(length-2), self.read(length-2))))

    def end_of_image(self):
       self.log(":END OF IMAGE")

    def decode(self):
        while self.i < self.n:
            ff = self.readb(1)[0]
            if ff != 0xFF: continue
            #assert ff == 0xFF, hex(ff)
            marker = self.readb(1)[0]
            self.j = self.i

            if marker == 0x00:
                continue
            elif marker == 0xD8:
                self.start_of_image()
            elif marker in [0xC0, 0xC1, 0xC2, 0xC3, 0xC9, 0xCA, 0xCB]:
                self.start_of_frame(marker)
            elif marker == 0xC4:
                self.huffman_table()
            elif marker == 0xDB:
                self.quantization_table()
            elif marker == 0xDA:
                self.start_of_scan()
            elif marker == 0xE0:
                self.application_zero()
            elif 0xE1 <= marker <= 0xEF:
                self.application_specific(marker)
            elif marker == 0xFE:
                self.comment()
            elif marker == 0xD9:
                self.end_of_image()
            else:
                assert True, marker
            print()

with open("god.jpg", "rb") as f:
    Decoder(f.read()).decode()
