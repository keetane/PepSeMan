import os

import uvicorn


if __name__ == "__main__":
    uvicorn.run(
        "app.api:app",
        host=os.getenv("PMMP_HOST", "127.0.0.1"),
        port=int(os.getenv("PMMP_PORT", "8000")),
        reload=False,
    )

