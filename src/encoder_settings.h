
/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2019 Pere Rafols Soler
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/


#ifndef RMSI_XBIN_ENCODER_SETTINGS_H
#define RMSI_XBIN_ENCODER_SETTINGS_H

#define IMG_STREAM_8bits //Comment this line out if you prefere to encode the imgStream as 16bits plus a mask

#ifdef IMG_STREAM_8bits
typedef  unsigned char imgstreamencoding_type;
#define ENCODER_RANGE 254.0 //255.0 is the maximum but I'm making it a bit below that to get some headroom
#define ENCODING_BITS 8
#define ENCODING_BIT_MASK 0xFF //imgStream encoding mask (in 8 bit encoding, no masking)
#else
typedef  unsigned short imgstreamencoding_type;
#define ENCODER_RANGE 65500.0 //65535.0 is the maximum but I'm making it a bit below that to get some headroom
#define ENCODING_BITS 16
#define ENCODING_BIT_MASK 0xFFC0 //imgStream encoding mask (only used for 16 bit encoding)
#endif


#endif