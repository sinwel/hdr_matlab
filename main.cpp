#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern void hdrprocess_sony_raw(unsigned short *src, unsigned short *dst, unsigned short *simage, int W, int H, int times, int noise_thred);
int main(int argc, char **argv)
{
	FILE			*fp;
	unsigned short	*src, *dst;
	const char		ifname[] = "G:\\Programs\\debug_files\\hdrTest1.raw";
	const char		ofname[] = "G:\\Programs\\debug_files\\hdrTestout.raw";
	int				W = 4032;
	int				H = 3024;

	fp = fopen(ifname, "rb+");
	src = (unsigned short*)malloc(W*H*sizeof(unsigned short));
	dst = (unsigned short*)malloc(W*H*sizeof(unsigned short));
	fread(src, 2, W*H, fp);
	fclose(fp);

	hdrprocess_sony_raw(src, dst, 0, W, H, 4, 0);

	fp = fopen(ofname, "wb+");
	fwrite(dst, 2, W*H, fp);
	fclose(fp);


	return 0;
}