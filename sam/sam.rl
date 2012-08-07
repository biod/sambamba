/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
%%{
    machine sam;
  
    newline = ('\n' | '\r' | '\r''\n');
    headervalue = [ -~]+;
    comment = "@CO\t" any*;
    
    user_tag = (lower alnum | upper lower) ;
    user_datafield = user_tag ':' headervalue ;

    format_version = digit+ '.' digit+ ;
    vn_field = "VN:" format_version ;
    sorting_order = "unknown" | "unsorted" | "queryname" | "coordinate" ;
    so_field = "SO:" sorting_order ;
    hd_datafield = vn_field | so_field ;
    hdline = "@HD" ('\t' (hd_datafield | user_datafield))+ ;

    sequence_name = [!-)+-<>-~] [!-~]* ;
    sn_field = "SN:" sequence_name ;
    sequence_length = [1-9] digit* ;
    ln_field = "LN:" sequence_length ;
    as_field = "AS:" headervalue ;
    md5 = [A-Fa-f0-9]+ ;
    m5_field = "M5:" md5 ;
    sp_field = "SP:" headervalue ;
    ur_field = "UR:" headervalue ;
    sq_datafield = sn_field | ln_field | as_field | m5_field | sp_field | ur_field ;
    sqline = "@SQ" ('\t' (sq_datafield | user_datafield))+ ;
 
    rgid_field = "ID:" headervalue ;
    cn_field = "CN:" headervalue ;
    ds_field = "DS:" headervalue ;
    dt_field = "DT:" headervalue ;
    fo_field = "FO:" ( '*' | [ACMGRSVTWYHKDBN]+ ) ;
    ks_field = "KS:" headervalue ;
    lb_field = "LB:" headervalue ;
    pg_field = "PG:" headervalue ;
    pi_field = "PI:" headervalue ;
    pl_field = "PL:" "CAPILLARY" | "LS454" | "ILLUMINA" | "SOLID" | 
                     "HELICOS" | "IONTORRENT" | "PACBIO" ;
    pu_field = "PU:" headervalue ;
    sm_field = "SM:" headervalue ;
    rg_datafield = rgid_field | cn_field | ds_field | dt_field | fo_field | ks_field |
                   lb_field | pg_field | pi_field | pl_field | pu_field | sm_field;
    rgline = "@RG" ('\t' (rg_datafield | user_datafield))+ ;

    pgid_field = "ID:" headervalue ;
    pn_field = "PN:" headervalue ;
    cl_field = "CLL" headervalue ;
    pp_field = "PP:" headervalue ;
    pgvn_field = "VN:" headervalue ;
    pg_datafield = pgid_field | pn_field | cl_field | pp_field | pgvn_field ;
    pgline = "@PG" ('\t' (pg_datafield | user_datafield))+ ;

    headerline = comment | hdline | sqline | rgline | pgline ;
    header = (hdline newline)? headerline (newline headerline)+ ;

    sign = [\-+] ;

    qname =  '*' | [!-?A-~]{1,255} ;
    flag = [0-9]{1,5} ;
    rname = '*' | [!-()+-<>-~] [!-~]* ;
    pos = [0-9]{1,9} ;
    mapq = [0-9]{1,3} ;
    cigar = '*' | (digit+ [MIDNSHPX=])+ ; # TODO: more validation
    rnext = '*' | '=' | [!-()+-<>-~][!-~]* ;
    pnext = digit{1,9} ;
    tlen = sign? digit{1,9} ;
    seq = '*' | [A-Za-z=.]+ ;
    qual = [!-~]+ ;

    mandatoryfields = qname '\t' 
                       flag '\t' 
                      rname '\t' 
                        pos '\t' 
                       mapq '\t'
                      cigar '\t'
                      rnext '\t'
                      pnext '\t'
                       tlen '\t'
                        seq '\t'
                       qual ;

    charvalue =  [!-~] ;
    integervalue = sign? digit+ ;
    floatvalue = sign? digit* '.'? digit+ ([eE] sign? digit+)? ; 
    stringvalue = [ !-~]+ ;
    hexstringvalue = xdigit+ ;
    integerarrayvalue = [cCsSiI] (',' integervalue)+ ;
    floatarrayvalue = [f] (',' floatvalue)+ ;
    arrayvalue = integerarrayvalue | floatarrayvalue ;

    tagvalue = ("A:" charvalue) | 
               ("i:" integervalue) | 
               ("f:" floatvalue) | 
               ("Z:" stringvalue) | 
               ("H:" hexstringvalue) |
               ("B:" arrayvalue) ;

    tag = alpha alnum ; # TODO: add predefined tags
    optionalfield = tag ':' tagvalue ;
    optionalfields = optionalfield ('\t' optionalfield)* ;

    alignment = mandatoryfields ('\t' optionalfields)? ;

    sam := header? newline (alignment (newline alignment)*)? ;

    write data; 
}%%
