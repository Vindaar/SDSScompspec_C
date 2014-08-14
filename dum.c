int main ()
{
int c, i;

// There are 256 ascii characters, from 0 to 255.
// A lot of the values just can't be displayed. For example, there is
// the newline character, a backspace character, a tab character, and
// probably the most entertaining one, the beep character.

for (c = 0; c < 256; )
{
// Adjust the stopping value for i accordingly.
for (i = 0; i < (80) / 7 && c < 256; i++, c++)
{
printf ("%3d %c ", c, c);
}
// Depending on stopping value, this may or may not be commented out.
printf ("\n");
}
return 0;
}
