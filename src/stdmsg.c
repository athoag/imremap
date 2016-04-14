//#ifndef lint
//static char vcid[] = "$Id: stdmsg.c,v 1.1.1.1 2008-12-17 23:23:09 marusa Exp $";
//#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "stdmsg.h"


/**************************************************************************/
/**************************************************************************/
//#define MSG_LINE_LENGTH 75
#define MSG_LINE_LENGTH 100
/**************************************************************************/
/**************************************************************************/
#define MSG_LENGTH (MSG_LINE_LENGTH - 15)

FILE *fmsg = NULL;

int teeprintf(const char *fmt, ...) 
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)))
#endif
     ;
void printmessage(const char *prompt, const char *fmt, va_list args)
#ifdef __GNUC__
     __attribute__ ((format (printf, 2, 0)))
#endif
     ;


int teeprintf(const char *fmt, ...)
{ 
  va_list args;
  int r;
  
  va_start(args, fmt);
  r = vfprintf(stderr, fmt, args);
  if (fmsg) {
    vfprintf(fmsg, fmt, args);
    fflush(fmsg);
  }
  va_end(args);
  return r;
}

/* msgseparator - Prints a separator line, similar to the one shown after
 * the caption.
 */
void msgseparator(void) {
  char line[MSG_LENGTH+1];
  int i;

  for (i=0; i<MSG_LENGTH; i++) line[i] = '-';
  line[MSG_LENGTH] = 0;
  teeprintf("-------------:-%s\n", line);
}
  

/* printclocktime - Prints on stderr the first 4 characters of prompt,
 * then the CPU time used so far by the program, and finally a column
 * (:).  The first execution of the procedure arranges the caption to
 * be printed.  For instance we have
 * `printclocktime("Info")':
 *
 *     CPU time | Action 
 * -------------:-------------------------------------
 * Info[0000.00]: _
 *
 * The cursor is left at position signed by _.
 */
void printclocktime(const char *prompt) {
  static clock_t CPUstart = 1, CPUold;
  static float acc;
  clock_t temp;
  
  if (CPUstart == 1) {
    teeprintf("\n    CPU time : Action \n");
    msgseparator();
    CPUstart = clock();
    CPUold = 0;
    acc = 0.0;
  }
  if (clock() != (clock_t)(-1)) {
    temp = clock() - CPUstart;
    if (temp < CPUold) {
#ifdef __GNUC__
      acc += 4294967296.0 / (float)CLOCKS_PER_SEC;
#else
      acc += pow(2.0, sizeof(temp)*8) / (float)CLOCKS_PER_SEC;
#endif
    }
    teeprintf("%-4s[%07.2f]: ", prompt,
	      (float)(temp)/(float)CLOCKS_PER_SEC + acc);
    CPUold = temp;
  } else teeprintf("%-4s[????.??]: ", prompt);
}

/* printmessage - Prints an output message.  The message is composed
 * of 4 characters used as prompt (see printclocktime), of the format
 * for the message itself, and of the list of argoments args.  The
 * message, if needed, is wrapped is several lines automatically; a
 * new line can in any case always be forced using the character `\n'.
 * A new line is automatically added at the end of the message, unless
 * the line ends with the ESC (\e) character (in this case, the \e is
 * removed and no new line is added).  This feature allows one to
 * break the process of writing a line into several stages (i.e., when
 * using \e at the end of a message, the next message just behaves as
 * if it was joined to the first one).  One can always reset the
 * status by calling printmessage with "\r" as format (this force a
 * reset at the next call of printmessage; hence it basically cancels
 * the \e request; see input for a use of this).  If the prompt is
 * NULL, no new prompt is shown.  Tabulations are also allowed: a
 * series of tabs (\t) on a line is first considered as an indication
 * for marks; if tabs are found again in following lines, they are
 * considered as goto-mark.  The first tabulation is also used as
 * paragraph indentation mark (hence, put a tab just after a newline
 * if you don't want lines to be indented).  Another fetaure is given
 * by the so-colled ghost-line.  If the text to print starts with a
 * line made only of spaces and tabs, this line is not printed but
 * rather used as `model' for tabs.  This is useful to set at the
 * beginning the tabulations to be used for the rest of the text.
 */
void printmessage(const char *prompt, const char *fmt, va_list args)
{ 
  char buffer[16384], *endbuffer, *par, *endpar, *word, *endword, *tabpos;
  int endnewline;
  static int line = 0, linetabs = 0, pos = 0, ghostmode = 0,
    ntabs = 0, tabs[16], curtab = 0,
    reset = 1;
  
  vsprintf(buffer, fmt, args);
  if (strcmp(buffer, "\r") == 0) {
    reset = 1;
    return;
  }
  if (buffer[strlen(buffer) - 1] == '\e') {
    endnewline = 0;
    buffer[strlen(buffer) - 1] = 0;
  } else endnewline = 1;
  if (reset) {
    if (prompt != NULL) printclocktime(prompt);
    else teeprintf("             : ");
  }
  endbuffer = strchr(buffer, (char)0);
  word = buffer;
  par = buffer;
  if (reset) line = linetabs = ntabs = 0;
  while (1) {
    if (reset) pos = curtab = ghostmode = 0;
    reset = 1;
    endpar = strchr(par, '\n');
    if (endpar == NULL) endpar = endbuffer;
    else if (word == buffer) {
      ghostmode = 1;
      for (word=par; word<endpar; word++) {
	if ((*word != ' ') && (*word != '\t')) {
	  ghostmode = 0;
	  break;
	}
      }
      word = buffer;
    }
    *endpar = (char)0;
    word = par;
    while (word < endpar) {
      endword = strchr(word, ' ');
      if (endword == NULL) endword = endpar;
      *endword = (char)0;
      if (*word == '\t') {
	word++;
	if ((ntabs == 0) || (line == linetabs)) {
	  if (ntabs < 16) tabs[ntabs++] = pos;
	  linetabs = line;
	} else {
	  if (curtab < ntabs) {
	    if (ghostmode) pos = tabs[curtab];
	    else for (; pos<tabs[curtab]; pos++) teeprintf(" ");
	    curtab++;
	  }
	}
      }
      tabpos = strchr(word, '\t');
      if (tabpos) *tabpos = (char)0;
      pos += strlen(word);
      if (pos > MSG_LENGTH) {
	line++;
	if (! ghostmode) teeprintf("\n             : ");
	ghostmode = 0;
	if (ntabs > 0) 
	  for (pos=0; pos<tabs[0]; pos++) teeprintf(" ");
	else pos = 0;
	pos += strlen(word);
      }
      if (! tabpos) {
	if (endword != endpar) {
	  pos++;
	  if (! ghostmode) teeprintf("%s ", word);
	} else {
	  if (! ghostmode) teeprintf("%s", word);
	}
	word = endword + 1;
      } else {
	if (! ghostmode) teeprintf("%s", word);
	*tabpos = '\t';
	word = tabpos;
	*endword = ' ';
      }
    }
    par = endpar + 1;
    if ((par < endbuffer) || ((endpar < endbuffer) && (endnewline == 0))) {
      if (! ghostmode) teeprintf("\n             : ");
    } else break;
    line++;
    ghostmode = 0;
  }
  if (endnewline) {
    reset = 1;
    if (pos > -1) teeprintf("\n");
  } else reset = 0;
  fflush(stderr);
}

/* message - Shows a message on stderr.  Parameters are much like
 * printf.
 */
void message(const char *fmt, ...)
{ 
  va_list args;
  
  va_start(args, fmt);
  printmessage("Info", fmt, args);
  errno = 0;
  va_end(args);
}

/* warning - Shows a warning on stderr.  Parameters are much like
 * printf.
 */
void warning(const char *fmt, ...)
{  
  va_list args;

  va_start(args, fmt);
  printmessage("Warn", fmt, args);
  va_end(args);
}

/* error - Shows an error in stderr and quit from the program (via
 * `exit').  The first parameter is the name of the calling procedure,
 * other parameters as in printf.
 */
void error(const char *proc, const char *fmt, ...)
{  
  va_list args; 
  char newfmt[1024];

  va_start(args, fmt);
  sprintf(newfmt, "ERROR in %s: %s", proc, fmt);
  printmessage("Err", newfmt, args);
  if (errno != 0) {
    teeprintf("             : %s", strerror(errno));
    teeprintf("\n");
    message("Quitting.");
    exit(errno);
  }
  exit(-1);
  va_end(args);
}

/* input - This function shows a prompt according to the provided
 * format and arguments, and wait for user's input.  The input is then
 * returned as a static string.
 */
char *input(const char *fmt, ...)
{  
  va_list args; 
  char newfmt[1024];
  static char buffer[1024];

  va_start(args, fmt);
  sprintf(newfmt, "%s\e", fmt);
  printmessage("User", newfmt, args);
  if (fgets(buffer, sizeof(buffer), stdin) == NULL) 
    return NULL;
  else {
    if (fmsg) {
      fprintf(fmsg, "%s", buffer);
      fflush(fmsg);
    }
    printmessage(NULL, "\r", args);
    return buffer;
  }
}
  
  
/* bp - Shows just the string "Breakpoint <n>", where <n> is the
 * integer passed as parameter.  To be used during debbuging.
 */
//void bp(int n)
//{ 
//  char newfmt[1024];
//  
//  sprintf(newfmt, "Breakpoint <%d>", n);
//  printmessage("Brk", newfmt, (va_list)NULL);
//}

/* progressbar - Shows a progress bar for a task completed at a fraction
 * f (f <- [0, 1]).  The procedure inizializes the various variables
 * automatically.  Remember to call the procedure with f >= 1.0 at the
 * end of the task.
 */
void progressbar(float f)
{
  static float f0 = -1.0;
  static struct timeval t0, t1;
  float eta;
  int i, j;
  char line1[MSG_LINE_LENGTH], line2[MSG_LINE_LENGTH], *c;

  if (! isatty(1)) {
    if (f >= 1.0) message("Done.");
    return;
  } else {
    if (f >= 1.0) {
      f0 = -1.0;
      strcpy(line1, "Done.");
      c = &line1[5];
      for (i=5; i<MSG_LENGTH; i++, c++) *c = ' ';
      *c = 0;
      message("%s",line1);  ///////<<<<<<<<============
      return;
    } else {
      if (f0 < 0) {
	gettimeofday(&t0, NULL);
	f0 = f;
      }
      if (gettimeofday(&t1, NULL) == 0) {
	eta = t1.tv_sec - t0.tv_sec + 
	  ((float)t1.tv_usec - (float)t0.tv_usec) / 1000000.0;
	if (eta < 1e6*(f - f0)) 
	  fprintf(stderr, "ETA [%07.01f]: ", (1 - f) * eta / (f - f0));
	else 
	  fprintf(stderr, "ETA [-----.-]: ");
      } else fprintf(stderr, "ETA [?????.?]: ");
      j = (MSG_LENGTH-5) / 2;
      c = line1;
      for (c=line1, i=0; i<j; i++, c++) *c = ' ';
      sprintf(c, "% 3d%%", (int)((f*100.0) + 0.5));
      i = strlen(line1);
      c = &line1[i];
      for (; i<MSG_LENGTH-1; i++, c++) *c = ' ';
      *(c++) = '|';
      *c = 0;
      j = (int)(f*(MSG_LENGTH-1) + 0.5);
      strncpy(line2, line1, j);
      line2[j] = 0;
      c = &line1[j];
      fprintf(stderr, "\x1B[7m%s\x1B[0m%s\r", line2, c);
      fflush(stderr);
    }
  }
}

/* showcounter - Shows a a bouncing star inside brackets [ * ]
 * toghether with some message.  Should be used to report the progress
 * of an operation for which no time estimation is available.  The
 * message can (and perhaps should...) change for each different call
 * of the procedure.  The procedure is initialized with a call with
 * NULL message; then the repeated calls are made; finally everything
 * is cleaned up with a call with empty ("") message.
 */
void showcounter(const char *fmt, ...)
{
  va_list args;
  char prompt[16], line[MSG_LENGTH+1];
  clock_t curCPU;
  static int p = 0, s = 1, lastlen = 0;
  static clock_t lastCPU = (clock_t)0;
  int i;
  
  if (! isatty(1)) return;
  if (fmt == NULL) {
    p = 0;
    s = 1;
    lastlen = 0;
    lastCPU = clock();
    return ;
  } else {
    va_start(args, fmt);
    strcpy(prompt, "Wait[       ]");
    prompt[5 + p] = '*';
    fprintf(stderr, "%s: ", prompt);
    vsprintf(line, fmt, args);
    line[MSG_LENGTH] = 0;
    fprintf(stderr, "%s", line);
    for (i=strlen(line); i<lastlen; i++) fprintf(stderr, " ");
    lastlen = strlen(line);
    fprintf(stderr, "\r");
    fflush(stderr);
    curCPU = clock();
    if ((curCPU == (clock_t)(-1)) ||
	(curCPU - lastCPU > CLOCKS_PER_SEC / 4)) {
      p += s;
      if (p == 7) {
	p = 5;
	s = -1;
      } else if (p == -1) {
	p = 1;
	s = +1;
      }
      lastCPU = curCPU;
    }
    va_end(args);
  }
}


/* parseErrorStr - Shows a string and, below it, a marker indicating the
 * position of an error.  Typically used to show parsing errors in reading
 * files.
 */
const char* parseErrorStr(const char *line, int pos)
{
  char buf1[MSG_LENGTH+1], buf2[MSG_LENGTH+1];
  static char buf3[MSG_LENGTH*2+10];
  long offset, len, i;

  len = strlen(line);
  if (line[len-1] == '\n') len--;
  if (len > MSG_LENGTH) {
    offset = pos - MSG_LENGTH/2;
    if (offset + MSG_LENGTH < len) offset = len - MSG_LENGTH;
    if (offset < 0) offset = 0;
  } else offset = 0;
  strncpy(buf1, line + offset, MSG_LENGTH);
  buf1[MSG_LENGTH] = 0;
  if (buf1[strlen(buf1)-1] == '\n') buf1[strlen(buf1)-1] = 0;
  if (offset != 0) {
    buf1[0] = '.';
    buf1[1] = '.';
    buf1[2] = '.';
  }
  if (offset + MSG_LENGTH < len) {
    buf1[MSG_LENGTH-1] = '.';
    buf1[MSG_LENGTH-2] = '.';
    buf1[MSG_LENGTH-3] = '.';
  }
  for (i=offset; i<pos; i++) buf2[i-offset] = ' ';
  buf2[pos-offset] = '^';
  buf2[pos-offset+1] = 0;
  sprintf(buf3, "%s\n%s", buf1, buf2);
  return buf3;

}

/* dateStr - Returns a string with the local time.
 */
const char* dateStr(void) 
{ 
  time_t t;
  
  time(&t);
  return asctime(localtime(&t));
}



