import os, sys

class Builder(object):

    """
    String Builder with rst specific helpers.

    Methods
    -------
    line: Add a single line of text
    title: Add a formatted title
    section: Add a formatted section
    classs: Add a class declaration
    attr: Add a attriute declaration
    func: Add a function declaration
    refclass: Add a class reference
    refattr: Add a attribute reference
    reffunc: Add a function reference
    rawline: Add a plaintext formatted line
    lines: Add multiple lines of text
    rawlines: Add multiple lines of plaintext

    """

    ct = ""
    def line(self, line = ""):
        self.ct = self.ct + line + "\n"
    def title(self, title):
        self.line(title)
        self.line("=" * len(title))
        self.line()
    def section(self, title):
        self.line(title)
        self.line("-" * len(title))
        self.line()
    def classs(self, name):
        self.line(".. js:class:: " + name)
        self.line()
    def attr(self, name):
        self.line(".. js:attribute:: " + name)
        self.line()
    def func(self, name):
        self.line(".. js:function:: " + name)
        self.line()
    def refclass(self, name):
        self.line(":js:class:`" + name + "`")
        self.line()
    def refattr(self, name):
        self.line(":js:attr:`" + name + "`")
        self.line()
    def reffunc(self, name):
        self.line(":js:func:`" + name + "`")
        self.line()
    def rawline(self, line):
        self.line("| " + line)
    def lines(self, lines):
        for line in lines:
            self.line(line)
    def rawlines(self, rawlines):
        for rawline in rawlines:
            self.rawline(rawline)
        self.line()
        
    def __str__(self):
        return self.ct
    

class Attributes(object):

    """
    Attribute Section Abstraction

    Parameters
    ----------
    title: Title of the section (usually Attributes)
    lines: Lines of comment text to parse


    Methods
    -------
    build: Write annotated comments to a builder

    """

    def __init__(self, title, lines):
        self.title = title
        self.attributes = dict()
        self.attrorder = []
        curr = ""
        for line in lines:
            if " : " in line:
                parts = line.split(" : ")
                name = parts[0]
                type = parts[1]

                self.attributes[name] = dict()
                self.attributes[name]["name"] = name
                self.attributes[name]["type"] = type
                self.attributes[name]["text"] = []
                self.attrorder.append(name)
                curr = name
            else:
                self.attributes[curr]["text"].append(line)

    def build(self, builder):
        builder.section(self.title)
        for name in self.attrorder:
            attr = self.attributes[name]
            builder.attr(attr["name"])
            builder.ct = builder.ct + "    "
            builder.refclass(attr["type"])
            builder.lines(attr["text"])
            builder.line()
        
class Methods(object):

    """
    Methods Section Abstraction

    Parameters
    ----------
    title: Title of the section
    lines: Lines of comment text to parse


    Methods
    -------
    build: Write annotated comments to a builder

    """

    def __init__(self, title, lines):
        self.title = title
        self.methods = dict()
        self.methorder = []
        curr = ""
        for line in lines:
            if line.endswith(")"):
                self.methods[line] = dict()
                self.methods[line]["ident"] = line
                self.methods[line]["text"] = []
                self.methorder.append(line)
                curr = line
            else:
                self.methods[curr]["text"].append(line)

    def build(self, builder):
        builder.section(self.title)
        for ident in self.methorder:
            md = self.methods[ident]
            builder.func(md["ident"])
            builder.lines(md["text"])
            builder.line()
    
class Plain(object):

    """
    Methods Section Abstraction

    Parameters
    ----------
    title: Title of the section
    lines: Lines of comment text to parse


    Methods
    -------
    build: Write comments to a builder


    """

    def __init__(self, title, lines):
        self.title = title
        self.lines = lines

    def build(self, builder):
        builder.section(self.title)
        builder.lines(self.lines)
    

class ClassPage(object):

    """
    ClassPage of documentation

    Parameters
    ----------
    lines: Lines of comment text to parse
    name: name of page
    sig: class signature
    filename: source filename
    line: line number in source file


    Methods
    -------
    write: Write rst documentation page into specified folder

    """

    def __init__(self, lines, name, sig, filename, line):
        self.name = name
        self.sig  = sig
        self.filename = filename
        self.line = line
        self.sections = []
        self.title = None

        area = "Description:"
        areas = dict()
        areas[area] = []
        areaorder = [area]
        prev_empty = True
        for line in lines:
            if self.title is None:
                if line is not "":
                    self.title = line.strip()
                continue
            if line.strip().endswith(":") and not line.strip().endswith("::") and not line.strip().startswith(":"):
                area = line.strip()
                if areas.get(area) is None:
                    areaorder.append(area)
                    areas[area] = []
                continue
            if line.strip() == "":
                if prev_empty:
                    continue #skip multiple newlines
                prev_empty = True
            prev_empty = False
            areas[area].append(line)
        
        for name in areaorder:
            lines = areas[name]
            name = name[0:len(name)-1]
            if "Attributes" in name:
                self.sections.append(Attributes(name, lines))
            elif "Methods" in name or "Event" in name or name == "Dispatch":
                self.sections.append(Methods(name, lines))
            else:
                self.sections.append(Plain(name, lines))
            
    def write(self, folder):
        print("Generating " + self.name)
        builder = Builder()

        builder.title(self.title)
        builder.classs(self.name + self.sig)
        builder.line("`{0}:{1} <https://github.com/STScI-SSB/pandeia/blob/master/ui/client/js/{0}#L{1}>`_".format(self.filename, self.line))
        builder.line()

        for sect in self.sections:
            sect.build(builder)

        builder.line()
        builder.section("Parents")
        t = ""
        sp = self.name.split(".")
        sp = sp[0:len(sp)-1]
        for i in sp:
            t = t + i
            builder.refclass(t)
            t = t + "."

        with open(os.path.join(folder, self.name + ".rst"), "w") as f:
            f.write("{0}".format(builder))
        return self.name

class FuncPage(object):

    """
    FuncPage of documentation

    Parameters
    ----------
    lines: Lines of comment text to parse
    name: name of page
    sig: class signature
    filename: source filename
    line: line number in source file


    Methods
    -------
    write: Write rst documentation page into specified folder

    """

    def __init__(self, lines, name, sig, filename, line):
        self.name = name
        self.sig  = sig
        self.filename = filename
        self.line = line
        self.sections = []

        area = "Description:"
        areas = dict()
        areas[area] = []
        areaorder = [area]
        prev_empty = True
        for line in lines:
            if line.strip().endswith(":") and not line.strip().endswith("::") and not line.strip().startswith(":"):
                area = line.strip()
                if areas.get(area) is None:
                    areaorder.append(area)
                    areas[area] = []
                continue
            if line.strip() == "":
                if prev_empty:
                    continue #skip multiple newlines
                prev_empty = True
            prev_empty = False
            areas[area].append(line)
        
        for name in areaorder:
            lines = areas[name]
            name = name[0:len(name)-1]
            self.sections.append(Plain(name, lines))
            
    def write(self, folder):
        print("Generating " + self.name)
        builder = Builder()

        builder.title(self.name)
        builder.func(self.name + self.sig)
        builder.line("`{0}:{1} <https://github.com/STScI-SSB/pandeia/blob/master/ui/client/js/{0}#L{1}>`_".format(self.filename, self.line))
        builder.line()

        for sect in self.sections:
            sect.build(builder)

        outname = self.filename.replace(".js", "") + "." + self.name

        with open(os.path.join(folder, outname + ".rst"), "w") as f:
            f.write("{0}".format(builder))
        return outname

class File(object):

    """
    Source File to be parsed

    Parameters
    ----------
    root: base directory of filename
    filename: name of file to parse


    Methods
    -------
    write: write all pages contained in source file to specified folder

    """

    def __init__(self, root, filename):
        self.classpages = []
        self.funcpages = []

        with open(os.path.join(root, filename)) as f:
            #classes
            comments = []
            linenum = 0
            is_comment = False
            is_ident   = False
            for line in f:
                linenum = linenum + 1
                if "/*" in line:
                    is_comment = True
                    if is_ident:
                        print("WARNING: comment not attached to class on line {0}".format(linenum))
                        is_ident = False
                        comments = []
                    continue
                if "*/" in line:
                    is_comment = False
                    is_ident   = True
                    continue
                if is_ident:
                    if line.startswith("pandeia"):
                        ident = line.split(" ")[0].strip()
                        sig   = line.split(" = ")[1].split("{")[0].strip().replace("function", "", 1)

                        self.classpages.append(ClassPage(comments, ident, sig, filename, linenum))

                        is_ident = False
                        comments = []
                    continue
                if is_comment:
                    line = line.replace(" *", "", 1)
                    line = line.replace("\n", "")
                    line = line.replace("\r", "")
                    comments.append(line)
                
        with open(root + filename) as f:
            #functions
            comments = []
            linenum = 0
            is_comment = False
            is_ident   = False
            for line in f:
                linenum = linenum + 1
                if line.strip().startswith("/*"):
                    is_comment = True
                elif line.strip().endswith("*/"):
                    is_comment = False
                elif is_comment:
                    line = line.replace(" *", "", 1)
                    line = line.replace("\n", "")
                    line = line.replace("\r", "")
                    comments.append(line)
                elif line.strip().startswith("//"):
                    line = line.replace("//", "", 1)
                    line = line.replace("\n", "")
                    line = line.replace("\r", "")
                    comments.append(line)
                elif line.startswith("function"):
                    line = line.replace("function", "", 1).split("{")[0].strip()
                    sp = line.split("(", 1)
                    ident = sp[0]
                    sig   = "(" + sp[1]
                    if len(comments) > 0:
                        self.funcpages.append(FuncPage(comments, ident, sig, filename, linenum))
                    else:
                        print("WARNING: Skipping function without comments {0}".format(ident + sig))
                else:
                    comments = []
                
    def write(self, folder):
        classnames = []
        funcnames = []
        for page in self.classpages:
            name = page.write(os.path.join(folder, "class"))
            classnames.append(name)

        for page in self.funcpages:
            name = page.write(os.path.join(folder, "func"))
            funcnames.append(name)

        return (classnames, funcnames)

def main():
    """

    Parse files in specified directory and write generated documentation
    into the specified out directory, along with a generated index page

    """

    if len(sys.argv) != 3:
        print("Expected usage: ")
        print("  " + sys.argv[0] + " path_to_js_folder output_folder")
        exit(-1)

    js_folder = sys.argv[1]
    out_folder = sys.argv[2]
    class_folder = os.path.join(out_folder, "class")
    func_folder = os.path.join(out_folder, "func")

    if not os.path.exists(class_folder):
        os.mkdirs(class_folder)
    if not os.path.exists(func_folder):
        os.mkdirs(func_folder)

    #Remove old docs
    for root, dirs, files in os.walk(out_folder):
        for filename in files:
            if filename.endswith(".rst"):
                os.remove(os.path.join(root, filename))

    parsed = []
    for root, dirs, files in os.walk(js_folder):
        for filename in files:
            if not filename.endswith(".js"):
                print("Skipping " + root + filename)
                continue
            print("Parsing " + root + filename)
            parsed.append(File(root, filename))


    classnames = []
    funcnames = []
    for f in parsed:
        cnames, fnames = f.write(out_folder)
        classnames = classnames + cnames
        funcnames = funcnames + fnames

    b = Builder()
    b.title("Pandeia Documentation")
    b.section("Classes")
    b.line()
    b.line(".. toctree::")
    b.line("  :maxdepth: 1")
    b.line()
    for name in sorted(classnames):
        b.line("  {0}<class/{0}>".format(name))
    b.line()

    b.section("Functions")
    b.line()
    b.line(".. toctree::")
    b.line("  :maxdepth: 1")
    b.line()
    for name in sorted(funcnames):
        b.line("  {0}<func/{0}>".format(name))

    with open(os.path.join(out_folder, "index.rst"), "w") as f:
        f.write("{0}".format(b))

main()
