import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public class Main extends JPanel {
    public static double gaussianModel(double x, double y, double variance) {
        return 1/(2*Math.PI*Math.pow(variance, 2)) * Math.exp((-x*x-y*y)/2/variance/variance);
    }

    public static double[][] generateWeightMatrix(int radius, double variance) {
        double[][] weights = new double[radius][radius];

        double sum = 0;
        for (int i = 0; i < radius; i++) {
            for (int j = 0; j < radius; j++) {
                weights[i][j] = gaussianModel(i-radius/2.0, j-radius/2.0, variance);
                sum += weights[i][j];
            }
        }

        for (int i = 0; i < radius; i++) {
            for (int j = 0; j < radius; j++) {
                weights[i][j] /= sum;
            }
        }

        return weights;
    }

    public static BufferedImage blur(BufferedImage img, double[][] kernel) {
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                double sum = 0;

                for (int k = 0; k < kernel.length; k++) {
                    for (int l = 0; l < kernel[k].length; l++) {

                        int x = i + k - (kernel.length / 2);
                        int y = j + l - (kernel.length / 2);

                        if (x >= img.getWidth() || y >= img.getHeight() || x < 0 || y < 0) {
                            continue;
                        }

                        Color sampled = new Color(img.getRGB(x, y));
                        int intensity = sampled.getRed();

                        sum += kernel[k][l] * intensity;
                    }
                }

                Color newColor = new Color(((int) sum), ((int) sum), ((int) sum));

                res.setRGB(i, j, newColor.getRGB());
            }
        }

        return res;
    }

    public static BufferedImage grayscale(BufferedImage orig) {
        BufferedImage res = new BufferedImage(orig.getWidth(), orig.getHeight(), BufferedImage.TYPE_INT_ARGB);
        for (int i = 0; i < orig.getWidth(); i++) {
            for (int j = 0; j < orig.getHeight(); j++) {
                int r = new Color(orig.getRGB(i, j)).getRed();
                int g = new Color(orig.getRGB(i, j)).getGreen();
                int b = new Color(orig.getRGB(i, j)).getBlue();

                res.setRGB(i, j, new Color((r+g+b)/3, (r+g+b)/3, (r+g+b)/3).getRGB());
            }
        }

        return res;
    }

    public static double[][] genSobelKernel(int k) {
        double[][] res = new double[k][k];

        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                double kc = (k-1)/2.0;
                res[i][j] = Math.abs(i-kc)/(i*i+j*j);
            }
        }

        return res;
    }

    private static double[][][] genGrads(BufferedImage img) {
        double[][] xkernel = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
        double[][] ykernel = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);
        double[][] gradx = new double[img.getWidth()][img.getHeight()], grady = new double[img.getWidth()][img.getHeight()];

        double maxx = -1, maxy = -1;

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                var bad = false;
                for(int a = 0; a < 3; a++) {
                    for(int b = 0; b < 3; b++) {
                        int x = i + a - 1;
                        int y = j + b - 1;

                        if (x >= img.getWidth() || y >= img.getHeight() || x < 0 || y < 0) {
                            gradx[i][j] = 0;
                            grady[i][j] = 0;
                            bad = true;
                        }

                        if (bad) break;

                        gradx[i][j] += new Color(img.getRGB(x, y)).getRed() * xkernel[a][b] / 256.0 / 4.0;
                        grady[i][j] += new Color(img.getRGB(x, y)).getRed() * ykernel[a][b] / 256.0 / 4.0;
                    }
                    if (bad) break;
                }

                maxx = Math.max(maxx, Math.abs(gradx[i][j]));
                maxy = Math.max(maxy, Math.abs(grady[i][j]));
            }
        }

        System.out.println(maxx + "     " + maxy);

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
//                gradx[i][j] *= 255/maxx;
//                gradx[i][j] = Math.min(gradx[i][j], 255);
//                grady[i][j] *= 255/maxy;
//                grady[i][j] = Math.min(grady[i][j], 255);

//                System.out.println("Before: " + gradx[i][j] + ", " + grady[i][j]);

//                gradx[i][j] = 255.0 / (1+Math.exp(0.7*(30-gradx[i][j])));
//                grady[i][j] = 255.0 / (1+Math.exp(0.7*(30-grady[i][j])));
//                halfway working

//                gradx[i][j] = 255.0 * (1 / (1+Math.exp(0.7*(30-255*gradx[i][j])))*2-1);
//                grady[i][j] = 255.0 * (1 / (1+Math.exp(0.7*(30-255*grady[i][j])))*2-1);

//                gradx[i][j] = 1.0 / (1+Math.exp(1/0.2*(30/255.0-gradx[i][j]))) * 2 - 1;
//                grady[i][j] = 1.0 / (1+Math.exp(1/0.2*(30/255.0-grady[i][j]))) * 2 - 1;

                gradx[i][j] = 1.0 / (1+Math.exp(-10*gradx[i][j])) * 2 - 1;
                grady[i][j] = 1.0 / (1+Math.exp(-10*grady[i][j])) * 2 - 1;

//                System.out.println("After: " + gradx[i][j] + ", " + grady[i][j]);
            }
        }

        return new double[][][] {gradx, grady};
    }

    public static BufferedImage sobelFilter(BufferedImage img) {
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);
        var grads = genGrads(img);
        var gradx = grads[0];
        var grady = grads[1];

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
//                System.out.println(gradx[i][j] + "    " + grady[i][j]);
                int mag = (int)(Math.sqrt((gradx[i][j]*gradx[i][j] + grady[i][j]*grady[i][j])/2) * 256);
                res.setRGB(i, j, new Color(mag,mag,mag).getRGB());
            }
        }

        return res;
    }

    public static BufferedImage sobelFilterWithNonMaxSuppress(BufferedImage img) {
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);
        var grads = genGrads(img);
        var gradx = grads[0];
        var grady = grads[1];

        double[][] newgradx = new double[img.getWidth()][img.getHeight()], newgrady = new double[img.getWidth()][img.getHeight()];

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                double dir = Math.atan2(grady[i][j], gradx[i][j]), dir2 = 0;
                for (double k = -Math.PI + Math.PI/8; k <= Math.PI; k += Math.PI/4) {
                    if (dir >= k - Math.PI/8 && dir < k + Math.PI/8) {
                        dir = k-Math.PI/8;
                        dir2 = k + Math.PI/8;
                        break;
                    }
                }

                int dx = (int) Math.signum(Math.cos(dir)), dy = (int) Math.signum(Math.sin(dir));
                double magnitude = Math.sqrt(gradx[i][j]*gradx[i][j] + grady[i][j]*grady[i][j]);
                boolean greatest = true;
                if (!(i + dx >= img.getWidth() || i + dx < 0 || j + dy >= img.getHeight() || j + dy < 0) &&
                        magnitude < Math.sqrt(gradx[i+dx][j+dy]*gradx[i+dx][j+dy] + grady[i+dx][j+dy]*grady[i+dx][j+dy]))
                    greatest = false;
                dx = (int) Math.signum(Math.cos(dir2));
                dy = (int) Math.signum(Math.sin(dir2));
                if (!(i + dx >= img.getWidth() || i + dx < 0 || j + dy >= img.getHeight() || j + dy < 0) &&
                        magnitude < Math.sqrt(gradx[i+dx][j+dy]*gradx[i+dx][j+dy] + grady[i+dx][j+dy]*grady[i+dx][j+dy]))
                    greatest = false;
//                if (!greatest || magnitude < 40/255.0) {
//                    newgradx[i][j] = 0;
//                    newgrady[i][j] = 0;
//                }
                if (!greatest) {
                    newgradx[i][j] = 0;
                    newgrady[i][j] = 0;
                }
                else {
                    newgradx[i][j] = gradx[i][j];
//                    newgradx[i][j] = 1;
                    newgrady[i][j] = grady[i][j];
//                    newgrady[i][j] = 1;

//                    newgradx[i][j] = 255 * ((1.0 / (1 + Math.exp(20-newgradx[i][j])))*2.0-1.0);
//                    newgrady[i][j] = 255 * ((1.0 / (1 + Math.exp(20-newgrady[i][j])))*2.0-1.0);
                }
            }
        }

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                int mag = (int)(Math.sqrt((newgradx[i][j]*newgradx[i][j] + newgrady[i][j]*newgrady[i][j])/2.0) * 255);
//                System.out.println(mag);
                res.setRGB(i, j, new Color(mag,mag,mag).getRGB());
            }
        }

        return res;
    }

    private final BufferedImage image;
    private final BufferedImage blurred;
    private final BufferedImage grayed;
    private final BufferedImage sobeled;
    private final BufferedImage nonmaxsobeled;

    public Main() {
        BufferedImage img = null;
        try {
            File imgFile = new File("src/img_10.png");
            System.out.println(imgFile.getAbsolutePath());
            img = ImageIO.read(imgFile);
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }
        assert img != null;
        this.image = img;

//        double radius = Math.sqrt(image.getWidth()*image.getHeight())/20;
        int radius = 5;
        var kernel = generateWeightMatrix(radius, Math.sqrt(radius));
        grayed = grayscale(img);
        blurred = blur(grayed, kernel);
        sobeled = sobelFilter(blurred);
        nonmaxsobeled = sobelFilterWithNonMaxSuppress(blurred);
//        nonmaxsobeled = null;

//        blurred = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_INT_ARGB);
//        for (int i = 0; i < blurred.getWidth(); i++) {
//            for (int j = 0; j < blurred.getHeight(); j++) {
//                blurred.setRGB(i, j, image.getRGB(i, j));
//            }
//        }

        try {
            ImageIO.write(blurred, "png", new File("src/blurred.png"));
            ImageIO.write(grayed, "png", new File("src/grayed.png"));
            ImageIO.write(sobeled, "png", new File("src/sobeled.png"));
            ImageIO.write(nonmaxsobeled, "png", new File("src/nonmaxsobeled.png"));
        } catch (IOException e) {
            e.printStackTrace();
        }

//        JFrame frame = new JFrame("image");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.setBounds(0, 0, 0, 0);
//        frame.setLayout(null);
//        Dimension prefSize = new Dimension(img.getWidth()*3, img.getHeight()*2);
//        frame.getContentPane().setPreferredSize(prefSize);
//        this.setBounds(0, 0, prefSize.width, prefSize.height);
//        System.out.println(img.getWidth() + " " + img.getHeight());
//        frame.add(this);
//        frame.pack();
//
//        frame.setVisible(true);
//        repaint();

        System.out.println("done");
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        setBackground(Color.BLACK);

        System.out.println("repainting");
        ImageIcon icon = new ImageIcon(image), icon2 = new ImageIcon(grayed), icon3 = new ImageIcon(blurred);
        ImageIcon icon4 = new ImageIcon(sobeled);
//        ImageIcon icon5 = new ImageIcon(nonmaxsobeled);
        g.drawImage(icon.getImage(), 0, 0, icon.getIconWidth(), icon.getIconHeight(), null);
        System.out.println("icon width: " + icon.getIconWidth());
        g.drawImage(icon2.getImage(), icon.getIconWidth(), 0, icon2.getIconWidth(), icon2.getIconHeight(), null);
        g.drawImage(icon3.getImage(), icon.getIconWidth() + icon2.getIconWidth(), 0, icon3.getIconWidth(), icon3.getIconHeight(), null);
        g.drawImage(icon4.getImage(), 0, icon.getIconHeight(), icon4.getIconWidth(), icon4.getIconHeight(), null);
//        g.drawImage(icon5.getImage(), icon4.getIconWidth(), icon2.getIconHeight(), icon5.getIconWidth(), icon5.getIconHeight(), null);
    }

    public static void main(String[] args) {
        var startTime = System.currentTimeMillis();
        new Main();
        var now = System.currentTimeMillis();
        System.out.println((now - startTime)/1000.0);
    }
}
