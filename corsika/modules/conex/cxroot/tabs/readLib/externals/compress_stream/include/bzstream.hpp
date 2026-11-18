// ============================================================================
// bzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2006  Frank Mattern 
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : bzstream.h
// Revision      : $Revision: 1.2 $
// Revision_date : $Date: 2009/05/28 11:36:30 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner, Frank Mattern
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================


// ============================================================================
// Update by Pierre Marijon <pierre@marijon.fr> in 2018 to become header only
// ============================================================================

#ifndef BZSTREAM_H
#define BZSTREAM_H 1


// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

#include <bzlib.h>

// ----------------------------------------------------------------------------
// Internal classes to implement bzstream. See below for user classes.
// ----------------------------------------------------------------------------

class bzstreambuf : public std::streambuf {
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for ibzstream at the end.

    BZFILE      *file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    char             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    bzstreambuf() : opened(0) {
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
    int is_open() { return opened; }
    bzstreambuf* open( const std::string& name, int open_mode);
    bzstreambuf* close();
    ~bzstreambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};

class bzstreambase : virtual public std::ios {
protected:
    bzstreambuf buf;
public:
    bzstreambase() { init(&buf); }
    bzstreambase( const std::string& name, int open_mode);
    ~bzstreambase();
    void open( const std::string& name, int open_mode);
    void close();
    bzstreambuf* rdbuf() { return &buf; }
};

// ----------------------------------------------------------------------------
// User classes. Use ibzstream and obzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the bz* 
// function interface of the zlib. Files are compatible with bzip compression.
// ----------------------------------------------------------------------------

/**
 * Input bzip stream class. 
 * Supports reading of bzip files. 
 * 
 * @author Deepak Bandyopadhyay, Lutz Kettner, Frank Mattern 
 */
class ibzstream : public bzstreambase, public std::istream {
public:
    ibzstream() : std::istream( &buf) {} 
    ibzstream( const std::string& name, int open_mode = std::ios::in)
        : bzstreambase( name, open_mode), std::istream( &buf) {}  
    bzstreambuf* rdbuf() { return bzstreambase::rdbuf(); }
    void open(  const std::string& name, int open_mode = std::ios::in) {
        bzstreambase::open( name, open_mode);
    }
    bool is_open() { return bzstreambase::rdbuf() ? bzstreambase::rdbuf()->is_open() : false; }
};

/**
 * Output bzip stream class. 
 * Supports writing of bzip files. 
 * 
 * @author Deepak Bandyopadhyay, Lutz Kettner, Frank Mattern 
 */
class obzstream : public bzstreambase, public std::ostream {
public:
    obzstream() : std::ostream( &buf) {}
    obzstream( const char* name, int mode = std::ios::out)
        : bzstreambase( name, mode), std::ostream( &buf) {}  
    bzstreambuf* rdbuf() { return bzstreambase::rdbuf(); }
    void open( const std::string& name, int open_mode = std::ios::out) {
        bzstreambase::open( name, open_mode);
    }
    bool is_open() { return bzstreambase::rdbuf() ? bzstreambase::rdbuf()->is_open() : false; }
};

// ----------------------------------------------------------------------------
// Internal classes to implement bzstream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class bzstreambuf:
// --------------------------------------

bzstreambuf* bzstreambuf::open( const std::string& name, int open_mode) {
    if ( is_open())
        return (bzstreambuf*)0;
    mode = open_mode;
    // no append nor read/write mode
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || ((mode & std::ios::in) && (mode & std::ios::out)))
        return (bzstreambuf*)0;
    char  fmode[10];
    char* fmodeptr = fmode;
    if ( mode & std::ios::in)
        *fmodeptr++ = 'r';
    else if ( mode & std::ios::out)
        *fmodeptr++ = 'w';
    *fmodeptr++ = 'b';
    *fmodeptr = '\0';
    file = 0; 
    file = BZ2_bzopen( name.c_str(), fmode);
    if (file == 0)
        return (bzstreambuf*)0;
    opened = 1;
    return this;
}

bzstreambuf * bzstreambuf::close() {
    if ( is_open()) {
        sync();
        opened = 0;
        BZ2_bzclose(file);
        return this;
    }
    return (bzstreambuf*)0;
}

int bzstreambuf::underflow() { // used for input buffer only
    if ( gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>( gptr());

    if ( ! (mode & std::ios::in) || ! opened)
        return EOF;
    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if ( n_putback > 4)
        n_putback = 4;
    std::memcpy( buffer + (4 - n_putback), gptr() - n_putback, n_putback);

    int num = 0; 
    num = BZ2_bzread( file, buffer+4, bufferSize-4);
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg( buffer + (4 - n_putback),   // beginning of putback area
          buffer + 4,                 // read position
          buffer + 4 + num);          // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());    
}

int bzstreambuf::flush_buffer() {
    // Separate the writing of the buffer from overflow() and
    // sync() operation.
    int w = pptr() - pbase();
    if ( BZ2_bzwrite( file, pbase(), w) != w)
        return EOF;
    pbump( -w);
    return w;
}

int bzstreambuf::overflow( int c) { // used for output buffer only
    if ( ! ( mode & std::ios::out) || ! opened)
        return EOF;
    if (c != EOF) {
        *pptr() = c;
        pbump(1);
    }
    if ( flush_buffer() == EOF)
        return EOF;
    return c;
}

int bzstreambuf::sync() {
    // Changed to use flush_buffer() instead of overflow( EOF)
    // which caused improper behavior with std::endl and flush(),
    // bug reported by Vincent Ricard.
    if ( pptr() && pptr() > pbase()) {
        if ( flush_buffer() == EOF)
            return -1;
    }
    return 0;
}

// --------------------------------------
// class bzstreambase:
// --------------------------------------

bzstreambase::bzstreambase( const std::string& name, int mode) {
    init( &buf);
    open( name, mode);
}

bzstreambase::~bzstreambase() {
    buf.close();
}

void bzstreambase::open( const std::string& name, int open_mode) {
    if ( ! buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
}

void bzstreambase::close() {
    if ( buf.is_open())
        if ( ! buf.close())
            clear( rdstate() | std::ios::badbit);
}

#endif // BZSTREAM_H
