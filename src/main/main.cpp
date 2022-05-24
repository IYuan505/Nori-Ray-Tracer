/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/parser.h>
#include <nori/scene.h>
#include <nori/camera.h>
#include <nori/block.h>
#include <nori/timer.h>
#include <nori/bitmap.h>
#include <nori/sampler.h>
#include <nori/integrator.h>
#include <nori/gui.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <filesystem/resolver.h>
#include <thread>

using namespace nori;

static int threadCount = -1;
static bool gui = true;
static bool iterativeMode = false;
static int startI = 0;

static void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block) {
    const Camera *camera = scene->getCamera();
    const Integrator *integrator = scene->getIntegrator();

    Point2i offset = block.getOffset();
    Vector2i size  = block.getSize();

    /* Clear the block contents */
    block.clear();

    /* For each pixel and pixel sample sample */
    for (int y=0; y<size.y(); ++y) {
        for (int x=0; x<size.x(); ++x) {
            for (uint32_t i=0; i<sampler->getSampleCount(); ++i) {
                Point2f pixelSample = Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
                Point2f apertureSample = sampler->next2D();

                /* Sample a ray from the camera */
                Ray3f ray;
                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);

                /* Compute the incident radiance */
                sampler->i = i;
                value *= integrator->Li(scene, sampler, ray);

                /* Store in the image block */
                block.put(pixelSample, value);
            }
        }
    }
}

static void render(Scene *scene, const std::string &filename) {
    const Camera *camera = scene->getCamera();
    Vector2i outputSize = camera->getOutputSize();
    scene->getIntegrator()->preprocess(scene);

    /* Create a block generator (i.e. a work scheduler) */
    BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE);

    /* Determine the filename of the output bitmap */
    std::string outputName = filename;
    size_t lastdot = outputName.find_last_of(".");
    if (lastdot != std::string::npos)
        outputName.erase(lastdot, std::string::npos);

    for (int i = startI; i < (int) std::ceil(outputSize.y() / (float) NORI_BLOCK_SIZE); ++i) {
        blockGenerator.setStart(Vector2i(0, i));
        int remaining_y = NORI_BLOCK_SIZE;
        if (i == (int) std::ceil(outputSize.y() / (float) NORI_BLOCK_SIZE) - 1
            && outputSize.y() % NORI_BLOCK_SIZE != 0) 
            remaining_y = outputSize.y() % NORI_BLOCK_SIZE;
        ImageBlock tmp_result(Vector2i(outputSize.x(), remaining_y), camera->getReconstructionFilter());
        tmp_result.clear();

        int startIdx = i * (int) std::ceil(outputSize.x() / (float) NORI_BLOCK_SIZE);
        int endIdx = (i + 1) * (int) std::ceil(outputSize.x() / (float) NORI_BLOCK_SIZE);

        /* Do the following in parallel and asynchronously */
        std::thread render_thread([&] {
            tbb::task_scheduler_init init(threadCount);

            cout << "Rendering .. ";
            cout.flush();
            Timer timer;

            tbb::blocked_range<int> range(startIdx, endIdx);

            auto map = [&](const tbb::blocked_range<int> &range) {
                /* Allocate memory for a small image block to be rendered
                   by the current thread */
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
                    camera->getReconstructionFilter());

                /* Create a clone of the sampler for the current thread */
                std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

                for (int i=range.begin(); i<range.end(); ++i) {
                    /* Request an image block from the block generator */
                    blockGenerator.next(block);

                    /* Inform the sampler about the block to be rendered */
                    sampler->prepare(block);

                    /* Render all contained pixels */
                    renderBlock(scene, sampler.get(), block);

                    /* The image block has been processed. Now add it to
                       the "big" block that represents the entire image */
                    tmp_result.put(block);
                }
            };

            /// Default: parallel rendering
            tbb::parallel_for(range, map);

            /// (equivalent to the following single-threaded call)
            // map(range);

            cout << "Part " << i << " done. (took " << timer.elapsedString() << ")" << endl;
        });

        /* Shut down the user interface */
        render_thread.join();

        Bitmap *tmp_bitmap = tmp_result.toBitmap();
        tmp_bitmap->saveEXR(outputName + "_tmp_" + std::to_string(i));
        delete tmp_bitmap;
        
        if (iterativeMode && i < (int) std::ceil(outputSize.y() / (float) NORI_BLOCK_SIZE) - 1) return;
    }

    /* Allocate memory for the entire output image and clear it */
    ImageBlock result(outputSize, camera->getReconstructionFilter());
    result.clear();
    /* Create a window that visualizes the partially rendered result */
    NoriScreen *screen = nullptr;
    if (gui) {
        nanogui::init();
        screen = new NoriScreen(result);
    }
    for (int i = 0; i < (int) std::ceil(outputSize.y() / (float) NORI_BLOCK_SIZE); ++i) {
        Bitmap bitmap(outputName + "_tmp_" + std::to_string(i) + ".exr");
        for (int y=0; y<bitmap.rows(); ++y)
            for (int x=0; x<bitmap.cols(); ++x)
                result.coeffRef(y + i * NORI_BLOCK_SIZE, x) << bitmap.coeff(y, x), 1;
    }

    /* Enter the application main loop */
    if (gui){
        nanogui::mainloop(50.f);
        delete screen;
        nanogui::shutdown();
    }
    std::unique_ptr<Bitmap> bitmap(result.toBitmap());
    bitmap->saveEXR(outputName);
    bitmap->savePNG(outputName);

    const Denoiser *denoiser = scene->getDenoiser();
    if (denoiser) {
        /* denoiser only applies to PNG (0-255) */
        Bitmap *denoised_bitmap = denoiser->denoise(result.toBitmap(), scene->getSampler());
        denoised_bitmap->saveDenoisedPNG(outputName + "_denoised");
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Syntax: " << argv[0] << " <scene.xml> [--no-gui] [--threads N]" <<  endl;
        return -1;
    }

    std::string sceneName = "";
    std::string exrName = "";

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);
        if (token == "-t" || token == "--threads") {
            if (i+1 >= argc) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }
            threadCount = atoi(argv[i+1]);
            i++;
            if (threadCount <= 0) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }

            continue;
        }
        else if (token == "--no-gui") {
            gui = false;
            continue;
        }
        else if (token == "--start") {
            startI = atoi(argv[i+1]);
            iterativeMode = true;
            i++;
            continue;
        }

        filesystem::path path(argv[i]);

        try {
            if (path.extension() == "xml") {
                sceneName = argv[i];

                /* Add the parent directory of the scene file to the
                   file resolver. That way, the XML file can reference
                   resources (OBJ files, textures) using relative paths */
                getFileResolver()->prepend(path.parent_path());
            } else if (path.extension() == "exr") {
                /* Alternatively, provide a basic OpenEXR image viewer */
                exrName = argv[i];
            } else {
                cerr << "Fatal error: unknown file \"" << argv[i]
                     << "\", expected an extension of type .xml or .exr" << endl;
            }
        } catch (const std::exception &e) {
            cerr << "Fatal error: " << e.what() << endl;
            return -1;
        }
    }

    if (exrName !="" && sceneName !="") {
        cerr << "Both .xml and .exr files were provided. Please only provide one of them." << endl;
        return -1;
    }
    else if (exrName == "" && sceneName == "") {
        cerr << "Please provide the path to a .xml (or .exr) file." << endl;
        return -1;
    }
    else if (exrName != "") {
        if (!gui) {
            cerr << "Flag --no-gui was set. Please remove it to display the EXR file." << endl;
            return -1;
        }
        try {
            Bitmap bitmap(exrName);
            ImageBlock block(Vector2i((int) bitmap.cols(), (int) bitmap.rows()), nullptr);
            block.fromBitmap(bitmap);
            nanogui::init();
            NoriScreen *screen = new NoriScreen(block);
            nanogui::mainloop(50.f);
            delete screen;
            nanogui::shutdown();
        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    }
    else { // sceneName != ""
        if (threadCount < 0) {
            threadCount = tbb::task_scheduler_init::automatic;
        }
        try {
            std::unique_ptr<NoriObject> root(loadFromXML(sceneName));
            /* When the XML root object is a scene, start rendering it .. */
            if (root->getClassType() == NoriObject::EScene)
                render(static_cast<Scene *>(root.get()), sceneName);
        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    }

    return 0;
}