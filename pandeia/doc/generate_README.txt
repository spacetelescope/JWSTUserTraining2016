Documentation For Our Documentation Generator (DFODG for short)

Our client code has plenty of good documentation sprinkled around functions and classes.

We now have a tool to lightly parse the documentation and code to generate Sphynx (rst) files
of formatted documentation.

Sections:
    A section is a way of describing a area of comments (sometimes specially formatted) which
    have a title.  Any line ending with ".*:" will be considered a section.

    Text within a section is formatted directly by sphynx and any commands therein will be interpreted.

    Cross linking is supported through the js domain:
    http://www.sphinx-doc.org/en/stable/domains.html#the-javascript-domain

    /*
     * Section:
     *   any text
     * between the first
     *       section and the
     * next is grouped
     *
     * Section2:
     *   some more text
     *   .. js:function: some_function
     */

Classes:
    A function is considered a class when the assignment follows the pattern:
    pandeia.some.path = function(some, optional, parameters)

    A comment block (/* */) above the class is considered the documentation for that class

    The comment block is formatted as folows:
    /*
     * Short Class Description/Title
     *  some more wordy
     *  description
     * 
     * Section A:
     *   section text
     *
     * Section B:
     *   some more
     *   section text
     */

    Special Cases:
        wordy description:
            Any misc text such as the woordy description is automatically
            placed into the Description section.

        custom sections:
            Attributes Section
                when section title is == Attributes
            
            Methods Section
                when section title is Methods, Dispatch, or contains Event
        

Functions:
    A line is considered a function when the assignment follows the pattern:
    //some comment or /* some comment */
    function somename(some,optional, parameters)

    A comment block (/* */) above the class is considered the documentation for that class
    A group of uninterupted single line comments will also be used if present instead

    The comment block is formatted as folows:
    /*
     *  some wordy
     *  description
     * 
     * Section A:
     *   section text
     *
     * Section B:
     *   some more
     *   section text
     */

    Special Cases:
        wordy description:
            Any misc text such as the woordy description is automatically
            placed into the Description section.

Usage:
    SPHYNX_DIR="some/path/to/contain/sphynx

    Setup Sphynx:
        cd $SPHYNX_DIR
        sphynx-quickstart

    Generate Docs:
        cd $PANDEIA_ROOT/src/pandeia/
        python doc/generate.py ./ui/client/js/ $SPHYNX_DIR/source/
        cd $SPHYNX_DIR
        make html

    View HTML:
        visit file:///$SPHYNX_DIR/build/html/index.html
