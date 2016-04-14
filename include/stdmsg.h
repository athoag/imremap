/* 	$Id: stdmsg.h,v 1.1.1.1 2008-12-17 23:23:09 marusa Exp $	 */

/* -*-mode: C; comment-column: 60; -*- */
#ifndef __STDMSG_H__
#define __STDMSG_H__

#include <stdio.h>

extern FILE *fmsg;

void msgseparator(void);
void message(const char *fmt, ...) 
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)))
#endif
     ;
void warning(const char *fmt, ...) 
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)))
#endif
     ;
void error(const char *proc, const char *fmt, ...) 
#ifdef __GNUC__
     __attribute__ ((noreturn, format (printf, 2, 3)))
#endif
     ;
void bp(int n);
void progressbar(float f);
void showcounter(const char *fmt, ...);
const char* parseErrorStr(const char *line, int pos);
const char* dateStr(void);


#endif
