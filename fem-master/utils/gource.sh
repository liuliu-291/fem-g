#! bin/zsh

gource --output-custom-log gmshfem_gource_log.txt --file-filter "^(contrib)|^(src.Contrib)|^(src.build)|^(src.Addon)"

sed -i '' -e "s/|ARoyer|/|Anthony Royer|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|Ismail BADIA|/|Ismail Badia|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|ismail badia|/|Ismail Badia|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|Mattesi|/|Vanessa Mattesi|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|XavierAdriaens|/|Xavier Adriaens|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|pmarchne|/|Philippe Marchner|/g" gmshfem_gource_log.txt
sed -i '' -e "s/|Marchner|/|Philippe Marchner|/g" gmshfem_gource_log.txt

gource --seconds-per-day 5e-2 --highlight-dirs --hide-filenames --font-size 40 --dir-font-size 20 --user-font-size 35 --highlight-users -1280x720 --logo id/GmshFEM-blackalpha-small.png --logo-offset 2284x1164 --date-format "%B, %Y" gmshfem_gource_log.txt -o - | ffmpeg -y -r 60 -f image2pipe -vcodec ppm -i - -vcodec libx265 -preset ultrafast -pix_fmt yuv420p -crf 26 -threads 0 -bf 0 gmshfem_history.mp4
